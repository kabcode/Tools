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
using ThresholdImageFilterType = itk::ThresholdImageFilter<ImageType>;

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout << "HUImageToAttenuationImage.exe InputImageFileName OutputImageFileName" << std::endl;
		return EXIT_FAILURE;
	}

	double WaterAttenuationAtEnergy = 0.1707; // at 100kev

	auto ImageFileReader = ImageFileReaderType::New();
	ImageFileReader->SetFileName(argv[1]);
		
	auto MultiplyWaterAttenuationImageFilter = MultiplyImageFilterType::New();
	MultiplyWaterAttenuationImageFilter->SetInput1(ImageFileReader->GetOutput());
	MultiplyWaterAttenuationImageFilter->SetConstant2(WaterAttenuationAtEnergy);

	auto MultiplyConstantImageFilter = MultiplyImageFilterType::New();
	MultiplyConstantImageFilter->SetInput1(MultiplyWaterAttenuationImageFilter->GetOutput());
	MultiplyConstantImageFilter->SetConstant2(1.0/1000);

	auto AddWaterAttenuationImageFilter = AddImageFilterType::New();
	AddWaterAttenuationImageFilter->SetInput1(MultiplyConstantImageFilter->GetOutput());
	AddWaterAttenuationImageFilter->SetConstant2(WaterAttenuationAtEnergy);

	auto ThresholdImageFilter = ThresholdImageFilterType::New();
	ThresholdImageFilter->SetInput(AddWaterAttenuationImageFilter->GetOutput());
	ThresholdImageFilter->SetLower(0);
	ThresholdImageFilter->SetOutsideValue(0);

	auto ImageFileWriter = ImageFileWriterType::New();
	ImageFileWriter->SetInput(ThresholdImageFilter->GetOutput());
	std::string OutputImageFilename(argv[2]);
	ImageFileWriter->SetFileName(OutputImageFilename);

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

