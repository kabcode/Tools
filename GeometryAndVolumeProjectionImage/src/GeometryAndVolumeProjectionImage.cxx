#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkExtractImageFilter.h"

#include "rtkThreeDCircularProjectionGeometryXMLFileReader.h"
#include "rtkCudaForwardProjectionImageFilter.h"

using PixelType = float;
static constexpr unsigned int Dim3 = 3;
static constexpr unsigned int Dim2 = 2;
using Image2DType = itk::CudaImage<PixelType, Dim2>;
using Image3DType = itk::CudaImage<PixelType, Dim3>;

using ImageReaderType = itk::ImageFileReader<Image3DType>;
using ImageWriterType = itk::ImageFileWriter<Image2DType>;
using ExtractImageFilterType = itk::ExtractImageFilter<Image3DType, Image2DType>;

using ThreeDProjectionGeometryReaderType = rtk::ThreeDCircularProjectionGeometryXMLFileReader;
using CudaForwardProjectionFilterType = rtk::CudaForwardProjectionImageFilter<Image3DType, Image3DType>;

void
PrintUsage(std::string ProgramName)
{
	std::cout << "USAGE: \n";
	std::cout << ProgramName << " IntensityImageFile ProjectionGeometryFile \n<DetectorSizeX DetectorSizeY DetectorSpacingX DetectorSpacingY>\n\n";
	std::cout << "OUTPUT" << " ProjectionImage with default size of 1000x1000 px and spacing 0.2 mm/px\n\n";
	std::cout << std::endl;
}

int main(int argc, char* argv[])
{
	/* input handling and error checking */
	if (argc < 3 || argc > 7)
	{
		PrintUsage(argv[0]);
		return EXIT_FAILURE;
	}

	/* Load intensity and label image */
	auto ImageReader = ImageReaderType::New();
	ImageReader->SetFileName(argv[1]);
	try
	{
		ImageReader->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	/* Load projection geometry */
	auto ProjectionGeometryReader = ThreeDProjectionGeometryReaderType::New();
	ProjectionGeometryReader->SetFilename(argv[2]);

	try
	{
		ProjectionGeometryReader->GenerateOutputInformation();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	/* Project volume onto image plane */
	auto Projection = Image3DType::New();
	Image3DType::IndexType Index;
	Index.Fill(0);
	Image3DType::SizeType Size;
	if(argc > 3)
	{
		Size[0] = atoi(argv[3]);
		Size[1] = atoi(argv[4]);
		Size[2] = 1;
	}
	else
	{
		Size[0] = 1000;
		Size[1] = 1000;
		Size[2] = 1;
	}

	Image3DType::SpacingType Spacing;
	if (argc > 3)
	{
		Spacing[0] = atoi(argv[5]);
		Spacing[1] = atoi(argv[6]);
		Spacing[2] = 1;
	}
	else
	{
		Spacing[0] = 0.2;
		Spacing[1] = 0.2;
		Spacing[2] = 0.2;
	}
	Projection->SetSpacing(Spacing);

	Image3DType::RegionType Region(Index, Size);
	Projection->SetRegions(Region);
	Projection->Allocate();
	Projection->FillBuffer(0);

	auto ForwardProjectionFilter = CudaForwardProjectionFilterType::New();
	ForwardProjectionFilter->SetGeometry(ProjectionGeometryReader->GetOutputObject());
	ForwardProjectionFilter->SetInput(Projection);
	ForwardProjectionFilter->SetInput(1, ImageReader->GetOutput());

	/* Extract 2D image from 3D image */
	auto ExtractionImageFilter = ExtractImageFilterType::New();
	ExtractionImageFilter->SetInput(ForwardProjectionFilter->GetOutput());
	ExtractionImageFilter->SetDirectionCollapseToIdentity();
	Image3DType::IndexType EIndex;
	EIndex.Fill(0);
	Image3DType::SizeType ESize;
	ESize = Projection->GetLargestPossibleRegion().GetSize();
	ESize[2] = 0;
	Image3DType::RegionType ERegion(EIndex, ESize);
	ExtractionImageFilter->SetExtractionRegion(ERegion);

	/* Write image*/
	auto ImageFileWriter = ImageWriterType::New();
	ImageFileWriter->SetInput(ExtractionImageFilter->GetOutput());
	ImageFileWriter->SetFileName("Projection.nrrd");

	try
	{
		ImageFileWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}