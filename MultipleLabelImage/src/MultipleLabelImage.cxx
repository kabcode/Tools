#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkThresholdImageFilter.h"

const unsigned int Dimension = 3;
using PixelType = float;
using ImageType = itk::Image<PixelType, Dimension>;

using ImageFileReaderType = itk::ImageFileReader<ImageType>;
using ImageFileWriterType = itk::ImageFileWriter<ImageType>;

using AddImageFilterType = itk::AddImageFilter<ImageType>;
using MultiplyImageFilterType = itk::MultiplyImageFilter<ImageType>;


int main(int argc, char *argv[])
{
	std::vector<ImageType::Pointer> ImageVector(argc - 1);
	
	for (auto i = 1; i < argc; ++i)
	{
		auto ImageFileReader = ImageFileReaderType::New();
		ImageFileReader->SetFileName(argv[i]);
		ImageVector[i - 1] = ImageFileReader->GetOutput();
		ImageVector[i - 1]->Update();
		ImageVector[i - 1]->DisconnectPipeline();
		ImageFileReader = nullptr;
	}

	std::vector<ImageType::Pointer> ImageVector2(argc - 1);

	for (auto i = 0; i < ImageVector.size(); ++i)
	{
		auto MultiplyConstantFilter = MultiplyImageFilterType::New();
		MultiplyConstantFilter->SetInput1(ImageVector[i]);
		MultiplyConstantFilter->SetConstant2(i+1);
		ImageVector2[i] = MultiplyConstantFilter->GetOutput();
		ImageVector2[i]->Update();
		ImageVector2[i]->DisconnectPipeline();
		MultiplyConstantFilter = nullptr;
	}
		
	auto AddImageFilter = AddImageFilterType::New();
	AddImageFilter->SetInput1(ImageVector2[0]);

	for (auto i = 1; i < ImageVector2.size(); ++i)
	{
		AddImageFilter->SetInput2(ImageVector2[i]);
		AddImageFilter->Update();
		AddImageFilter->SetInput1(AddImageFilter->GetOutput());
	}

	auto ThresholdImageFilter = itk::ThresholdImageFilter<ImageType>::New();
	ThresholdImageFilter->SetInput(AddImageFilter->GetOutput());
	ThresholdImageFilter->SetUpper(argc);
	ThresholdImageFilter->SetOutsideValue(0);

	auto ImageFileWriter = ImageFileWriterType::New();
	ImageFileWriter->SetInput(ThresholdImageFilter->GetOutput());
	ImageFileWriter->SetFileName("LabelImage.nrrd");

	try
	{
		ImageFileWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
	}

	return EXIT_SUCCESS;
}