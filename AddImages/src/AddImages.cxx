#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAddImageFilter.h"
#include "itkResampleImageFilter.h"

const unsigned int Dimension = 3;
using PixelType = float;
using ImageType = itk::Image<PixelType, Dimension>;


using ImageFileReaderType = itk::ImageFileReader<ImageType>;
using ImageFileWriterType = itk::ImageFileWriter<ImageType>;

auto WriteImage(const ImageType::Pointer, const std::string&) -> void;

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout << argv[0] << " InputFile0 InputFile2 OutputFile" << std::endl;
		return EXIT_FAILURE;
	}

    auto ImageFileReader0 = ImageFileReaderType::New();
    ImageFileReader0->SetFileName(argv[1]);

	auto ImageFileReader1 = ImageFileReaderType::New();
	ImageFileReader1->SetFileName(argv[2]);

	auto AddImageFilter = itk::AddImageFilter<ImageType>::New();
	AddImageFilter->SetInput1(ImageFileReader0->GetOutput());
	AddImageFilter->SetInput2(ImageFileReader1->GetOutput());

	try
	{
		AddImageFilter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}
  

    WriteImage(AddImageFilter->GetOutput(), argv[3]);
	
	return EXIT_SUCCESS;
}

void
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