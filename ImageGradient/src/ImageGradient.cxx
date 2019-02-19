#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

const unsigned int Dimension = 3;
using PixelType = float;
using ImageType = itk::Image<PixelType, Dimension>;
using VectorPixelType = itk::CovariantVector<float>;
using VectorImageType = itk::Image<VectorPixelType, Dimension>;
using VectorImagePointer = VectorImageType::Pointer;

using ImageFileReaderType = itk::ImageFileReader<ImageType>;
using ImageFileWriterType = itk::ImageFileWriter<VectorImageType>;

using GradientImageFilterType = itk::GradientImageFilter<ImageType>;
using GradientMagnImageFilterType = itk::GradientMagnitudeImageFilter<ImageType, ImageType>;

void WriteGradientImage(VectorImagePointer gradimage, std::string outfilename);

enum METHOD
{
	GradientImageFilter

};

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout << argv[0] << " InputImageFileName -o OutputFilename <-m []> <-magn>" << std::endl;
		return EXIT_FAILURE;
	}

	std::string OutputFilename("");
	auto GenerateMagnitudeImage = false;
	auto GradientFilterMethod = METHOD::GradientImageFilter; //default method
	for (auto i = 1; i < argc; ++i)
	{
		if (std::string(argv[i]) == "-o")
		{
			++i;
			OutputFilename = argv[i];
		}
		if (std::string(argv[i]) == "-m")
		{
			++i;
			if (argv[i] == "gfm")
				GradientFilterMethod = GradientImageFilter;
		}
		if (std::string(argv[i]) == "-magn")
		{
			++i;
			GenerateMagnitudeImage = true;
		}
	}

	auto ImageFileReader = ImageFileReaderType::New();
	ImageFileReader->SetFileName(argv[1]);
	ImageFileReader->Update();

	if(GradientFilterMethod == GradientImageFilter)
	{
		auto GradientImageFilter = GradientImageFilterType::New();
		GradientImageFilter->SetInput(ImageFileReader->GetOutput());
		GradientImageFilter->Update();
		WriteGradientImage(GradientImageFilter->GetOutput(), OutputFilename);
	}

	if(GenerateMagnitudeImage)
	{
		auto GaussianSmoothingFilter = itk::RecursiveGaussianImageFilter<ImageType>::New();
		GaussianSmoothingFilter->SetInput(ImageFileReader->GetOutput());
		GaussianSmoothingFilter->SetSigma(1);

		auto GradientMagnImageFilter = GradientMagnImageFilterType::New();
		GradientMagnImageFilter->SetInput(GaussianSmoothingFilter->GetOutput());
		
		auto ImageFileWriter = itk::ImageFileWriter<ImageType>::New();
		ImageFileWriter->SetFileName("Magn" + OutputFilename);
		ImageFileWriter->SetInput(GradientMagnImageFilter->GetOutput());
		ImageFileWriter->Update();
	}
	
	
	return EXIT_SUCCESS;
}

void
WriteGradientImage(VectorImagePointer gradimage, std::string outfilename)
{
	auto ImageFileWriter = ImageFileWriterType::New();
	ImageFileWriter->SetInput(gradimage);
	std::string OutputImageFilename(outfilename);
	ImageFileWriter->SetFileName(OutputImageFilename);

	try
	{
		ImageFileWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
	}
}