#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMaskImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkRegionOfInterestImageFilter.h"

using PixelType = double;
using LabelType = unsigned int;
static constexpr unsigned int Dim = 3;
using ImageType = itk::Image<PixelType, Dim>;
using LabelImageType = itk::Image< LabelType, Dim>;

using ImageReaderType = itk::ImageFileReader<ImageType>;
using LabelReaderType = itk::ImageFileReader<LabelImageType>;
using ImageWriterType = itk::ImageFileWriter<ImageType>;

using Labeltype = int;
using ShapeLabelObjectType = itk::ShapeLabelObject<Labeltype, LabelImageType::ImageDimension>;
using LabelMapType = itk::LabelMap<ShapeLabelObjectType>;
using LabelMapTypePointer = LabelMapType*;
using LabelImageToShapeLabelMapFilterType = itk::LabelImageToShapeLabelMapFilter<LabelImageType>;

using RegionOfInterestFilterType = itk::RegionOfInterestImageFilter<ImageType, ImageType>;

void
CreateMaskImage(ImageType::Pointer, LabelImageType::Pointer, ImageType::Pointer, LabelType);

void
PrintUsage(std::string ProgramName)
{
	std::cout << "USAGE: \n";
	std::cout << ProgramName << " IntensityImage LabelImage\n\n";
	std::cout << "OUTPUT" << " Label_X.nrrd with X = {0....n} label in LabelImage\n\n";
	std::cout << std::endl;
}

int main(int argc, char* argv[])
{
	/* input handling and error checking */
	if (argc < 3)
	{
		PrintUsage(argv[0]);
		return EXIT_FAILURE;
	}

	/* Load intensity and label image */
	auto IntensityImageReader = ImageReaderType::New();
	IntensityImageReader->SetFileName(argv[1]);
	try
	{
		IntensityImageReader->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	auto LabelImageReader = LabelReaderType::New();
	LabelImageReader->SetFileName(argv[2]);

	/* Compute shape label map */
	auto LabelImageToShapeLabelMapFilter = LabelImageToShapeLabelMapFilterType::New();
	LabelImageToShapeLabelMapFilter->SetInput(LabelImageReader->GetOutput());

	try
	{
		LabelImageToShapeLabelMapFilter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	auto LabelMapPointer = LabelImageToShapeLabelMapFilter->GetOutput();
	std::cout << "Size: " << LabelMapPointer->GetLargestPossibleRegion().GetSize() << std::endl;

	/* Divide label image into single label images */
	auto numberOfLabel = LabelMapPointer->GetNumberOfLabelObjects();
	if (numberOfLabel == 0)
	{
		std::cout << "No label found." << std::endl;
		return EXIT_FAILURE;
	}
	
	for (auto i = 0; i < numberOfLabel; ++i)
	{
		auto LabelObject = LabelMapPointer->GetNthLabelObject(i);
		
		/* Mask the intensity image */
		auto MaskedImage = ImageType::New();
		CreateMaskImage(MaskedImage, LabelImageReader->GetOutput(), IntensityImageReader->GetOutput(), LabelObject->GetLabel());

		/* Crop masked intensity images to the bounding box */
		auto RegionOfInterestFilter = RegionOfInterestFilterType::New();
		RegionOfInterestFilter->SetInput(MaskedImage);
		RegionOfInterestFilter->SetRegionOfInterest(LabelObject->GetBoundingBox());

		/* Write image */
		std::string FileName("Label_");
		FileName.append(std::to_string(i));
		FileName.append(".nrrd");
		auto ImageWriter = ImageWriterType::New();
		ImageWriter->SetInput(RegionOfInterestFilter->GetOutput());
		ImageWriter->SetFileName(FileName);

		try
		{
			ImageWriter->Update();
		}
		catch (itk::ExceptionObject &EO)
		{
			EO.Print(std::cout);
			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;
}

void
CreateMaskImage(
	ImageType::Pointer maskimage,
	LabelImageType::Pointer labelimage,
	ImageType::Pointer intensityimage,
	LabelType label)
{
	auto region = intensityimage->GetLargestPossibleRegion();
	maskimage->SetRegions(region);
	maskimage->Allocate();
	maskimage->FillBuffer(0);

	maskimage->SetSpacing(intensityimage->GetSpacing());
	maskimage->SetDirection(intensityimage->GetDirection());
	maskimage->SetOrigin(intensityimage->GetOrigin());

	itk::ImageRegionIterator<ImageType> MaskIterator(maskimage, maskimage->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<ImageType> IntensityIterator(intensityimage, intensityimage->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<LabelImageType> LabelIterator(labelimage, labelimage->GetLargestPossibleRegion());

	/* Copy the intensity from one image to another */
	while (!MaskIterator.IsAtEnd())
	{
		if(LabelIterator.Get() == label)
		{
			MaskIterator.Set(IntensityIterator.Get());
		}
		++MaskIterator;
		++IntensityIterator;
		++LabelIterator;
	}
}