// STL includes
#include <filesystem>

// ITK includes
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkEuler3DTransform.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkResampleImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkExtractImageFilter.h"


using PixelType = float;
const unsigned int Dim3D = 3;
const unsigned int Dim2D = 2;
using ImageType = itk::Image<PixelType, Dim3D>;

using ImageReaderType = itk::ImageFileReader<ImageType>;
using ImageWriterType = itk::ImageFileWriter<ImageType>;
using TransformReaderType = itk::TransformFileReader;
using TransformWriterType = itk::TransformFileWriter;
using OutputImageWriterType = itk::ImageFileWriter<ImageType>;
using CenteredEulerTransformType = itk::CenteredEuler3DTransform<double>;
using EulerTransformPointer = CenteredEulerTransformType::Pointer;
using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
using ChangeImageInformationFilterType = itk::ChangeInformationImageFilter<ImageType>;

void
PrintUsage(char * programmname)
{
	std::cout << "\nInformation to " << programmname << std::endl;
	std::cout << "Usage: " << std::endl;
	std::cout << programmname << " VolumeFile <optional parameters>" << std::endl;
	std::cout << "\noptional parameters:\n" << std::endl;
	std::cout << "-o OutputFilename" << std::endl;
	std::cout << std::endl;
}

int main(int argc, char* argv[])
{
	// check input arguments
	if (argc < 3)
	{
		PrintUsage(argv[0]);
		return EXIT_FAILURE;
	}

	auto invers = false;

	// input handling
	std::string OutputFilename("TransformedObject.nrrd");
	for (auto i = 0; i < argc; ++i)
	{
		if( std::string(argv[i]) == "-o")
		{
			++i;
			OutputFilename = argv[i];
		}
		if (std::string(argv[i]) == "-inv")
		{
			invers = true;
		}
	}

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

	
	auto T = CenteredEulerTransformType::New();
	if(!invers)
	{
		T->SetRotation(itk::Math::pi / 2, 0, 0);
	}
	else
	{
		T->SetRotation(-itk::Math::pi / 2, 0, 0);
	}
		
	auto ChangeImageInformationFilter = ChangeImageInformationFilterType::New();
	ChangeImageInformationFilter->SetInput(ImageReader->GetOutput());
	ChangeImageInformationFilter->ChangeDirectionOn();
	ChangeImageInformationFilter->SetOutputDirection(T->GetMatrix() * ImageReader->GetOutput()->GetDirection());
	ChangeImageInformationFilter->ChangeOriginOn();
	ChangeImageInformationFilter->SetOutputOrigin(T->GetMatrix() * ImageReader->GetOutput()->GetOrigin());

	auto ImageWriter = ImageWriterType::New();
	ImageWriter->SetInput(ChangeImageInformationFilter->GetOutput());
	ImageWriter->SetFileName(OutputFilename);
	
	try
	{
		ImageWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
	
}
