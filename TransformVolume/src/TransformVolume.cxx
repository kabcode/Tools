// STL includes
#include <filesystem>

// ITK includes
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkEuler3DTransform.h"
#include "itkResampleImageFilter.h"
#include "itkExtractImageFilter.h"


using PixelType = float;
const unsigned int Dim3D = 3;
const unsigned int Dim2D = 2;
using InputImageType = itk::Image<PixelType, Dim3D>;
using OutputImageType = itk::Image<PixelType, Dim3D>;

using ImageReaderType = itk::ImageFileReader<InputImageType>;
using ImageWriterType = itk::ImageFileWriter<OutputImageType>;
using TransformReaderType = itk::TransformFileReader;
using TransformWriterType = itk::TransformFileWriter;
using OutputImageWriterType = itk::ImageFileWriter<OutputImageType>;
using EulerTransformType = itk::Euler3DTransform<double>;
using EulerTransformPointer = EulerTransformType::Pointer;
using ResampleImageFilterType = itk::ResampleImageFilter<InputImageType, OutputImageType>;

void
PrintUsage(char * programmname)
{
	std::cout << "\nInformation to " << programmname << std::endl;
	std::cout << "Usage: " << std::endl;
	std::cout << programmname << " VolumeFile TransformFile" << std::endl;
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

	auto TransformReader = TransformReaderType::New();
	TransformReader->SetFileName(argv[2]);
	try
	{
		TransformReader->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	auto TransformList = TransformReader->GetTransformList();
	std::cout << "Number of transforms = " << TransformList->size() << std::endl;
	auto Transformation = TransformList->begin()->GetPointer();
	Transformation->Print(std::cout);
	auto T = EulerTransformType::New();
	T->SetParameters(Transformation->GetParameters());
	T->SetFixedParameters(Transformation->GetFixedParameters());
	
	auto ResampleFilter = ResampleImageFilterType::New();
	ResampleFilter->SetInput(ImageReader->GetOutput());
	ResampleFilter->SetTransform(T);
	ResampleFilter->SetDefaultPixelValue(1);

	auto SampleVolume = OutputImageType::New();
	OutputImageType::IndexType idx;
	idx.Fill(0);
	auto sz = ImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
	OutputImageType::RegionType reg(idx, sz);
	SampleVolume->SetRegions(reg);
	SampleVolume->Allocate();
	auto or = ImageReader->GetOutput()->GetOrigin();
	or -= T->GetTranslation();
	EulerTransformType::MatrixType M = T->GetMatrix().GetInverse();
	auto or1 = M * or ;
	SampleVolume->SetOrigin(or1);
	auto dir = ImageReader->GetOutput()->GetDirection();
	dir *= T->GetMatrix().GetInverse();
	SampleVolume->SetDirection(dir);
	auto spc = ImageReader->GetOutput()->GetSpacing();
	SampleVolume->SetSpacing(spc);

	ResampleFilter->SetReferenceImage(SampleVolume);
	ResampleFilter->UseReferenceImageOn();

	try
	{
		ResampleFilter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	auto ImageWriter = ImageWriterType::New();
	ImageWriter->SetInput(ResampleFilter->GetOutput());
	ImageWriter->SetFileName("TansformedObject.nrrd");

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
