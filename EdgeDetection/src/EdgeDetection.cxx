#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"

const unsigned int Dimension = 3;
using PixelType = float;
using ImageType = itk::Image<PixelType, Dimension>;

using ImageFileReaderType = itk::ImageFileReader<ImageType>;
using ImageFileWriterType = itk::ImageFileWriter<ImageType>;

using CannyEdgeFilterType = itk::CannyEdgeDetectionImageFilter<ImageType, ImageType>;


int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout << argv[0] << " InputImageFileName -o OutputFilename " << std::endl;
		return EXIT_FAILURE;
	}

	std::string OutputFilename("");
	auto GenerateMagnitudeImage = false;

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
			if (argv[i] == "gfm");
			
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

	auto CannyImageFilter = CannyEdgeFilterType::New();
	CannyImageFilter->SetInput(ImageFileReader->GetOutput());
	CannyImageFilter->SetLowerThreshold(10000);

	auto ImageFileWriter = ImageFileWriterType::New();
	ImageFileWriter->SetInput(CannyImageFilter->GetOutput());
	ImageFileWriter->SetFileName(OutputFilename);
	ImageFileWriter->Update();


	return EXIT_SUCCESS;
}