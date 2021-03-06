#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"

const unsigned int Dimension = 3;
using PixelType = float;
using ImageType = itk::Image<PixelType, Dimension>;


using ImageFileReaderType = itk::ImageFileReader<ImageType>;
using ImageFileWriterType = itk::ImageFileWriter<ImageType>;

auto ReadImage(const std::string&)->ImageType::Pointer;
auto WriteImage(const ImageType::Pointer, const std::string&) -> void;

int main(int argc, char *argv[])
{
	if(argc < 6)
	{
		std::cout << argv[0] << " InputImageFileName -o OutputFilename -spacing x y z" << std::endl;
		return EXIT_FAILURE;
	}

    auto ImageFileReader = ImageFileReaderType::New();
    ImageFileReader->SetFileName(argv[1]);
    auto ResamplingImage = ImageFileReader->GetOutput();
    try
    {
        ResamplingImage->Update();
    }
    catch (itk::ExceptionObject &EO)
    {
        EO.Print(std::cout);
        return EXIT_FAILURE;
    }
 

    std::string OutputFilename{ "" };
    auto outputSpacing = ImageFileReader->GetOutput()->GetSpacing();
	for (auto i = 1; i < argc; ++i)
	{
		if (std::string(argv[i]) == "-o")
		{
			++i;
			OutputFilename = argv[i];
		}
        if (std::string(argv[i]) == "-spacing")
        {
            ++i;
            outputSpacing[0] = CheckSpacing(std::strtof(argv[i],nullptr), outputSpacing[0]);
            ++i;
            outputSpacing[1] = CheckSpacing(std::strtof(argv[i], nullptr), outputSpacing[1]);
            ++i;
            outputSpacing[2] = CheckSpacing(std::strtof(argv[i], nullptr), outputSpacing[2]);
        }
	}

    auto ResampleFilter = itk::ResampleImageFilter<ImageType, ImageType>::New();
    ResampleFilter->SetInput(ResamplingImage);
  
    ResampleFilter->SetOutputDirection(ImageFileReader->GetOutput()->GetDirection());
    ResampleFilter->SetOutputOrigin(ImageFileReader->GetOutput()->GetOrigin());
    const auto OutputSize = ResizeDimension(ResamplingImage->GetLargestPossibleRegion().GetSize(), ResamplingImage->GetSpacing(), outputSpacing);
    ResampleFilter->SetSize(OutputSize);
    ResampleFilter->SetOutputSpacing(outputSpacing);
    ResampleFilter->Update();

    WriteImage(ResampleFilter->GetOutput(), OutputFilename);
	
	return EXIT_SUCCESS;
}

auto
ReadImage(const std::string&)->ImageType::Pointer
{
    
}

auto
WriteImage(const ImageType::Pointer image, const std::string& outfilename)
{
	auto ImageFileWriter = ImageFileWriterType::New();
	ImageFileWriter->SetInput(image);
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