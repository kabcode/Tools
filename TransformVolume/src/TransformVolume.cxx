// STL includes
#include <filesystem>

// ITK includes
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkEuler3DTransform.h"
#include "itkBSplineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkExtractImageFilter.h"


using PixelType = float;
const unsigned int Dim3D = 3;
const unsigned int Dim2D = 2;
using InputImageType = itk::Image<PixelType, Dim3D>;
using OutputImageType = itk::Image<PixelType, Dim3D>;

using ImageReaderType = itk::ImageFileReader<InputImageType>;
using ImageWriterType = itk::ImageFileWriter<OutputImageType>;
using TransformReaderType = itk::TransformFileReaderTemplate<double>;
using TransformWriterType = itk::TransformFileWriterTemplate<double>;
using OutputImageWriterType = itk::ImageFileWriter<OutputImageType>;
using EulerTransformType = itk::Euler3DTransform<double>;
using EulerTransformPointer = EulerTransformType::Pointer;
using BSplineTransformType = itk::BSplineTransform<double>;
using BSplineTransformPointer = BSplineTransformType::Pointer;
using ResampleImageFilterType = itk::ResampleImageFilter<InputImageType, OutputImageType>;

double inline deg2rad(double deg) { return deg / 180 * (itk::Math::pi / 2);  }

void
PrintUsage(char * programmname)
{
	std::cout << "\nInformation to " << programmname << "\n";
	std::cout << "Usage: " << "\n";
	std::cout << programmname << " VolumeFile -o OutputFile <optional parameters>" << "\n";
	std::cout << "\noptional parameters:\n" << "\n";
	std::cout << "-t x y z \t\t Translation in x,y,z direction in mm" << "\n";
	std::cout << "-r x y z \t\t Rotation around x,y,z axis in degree(order: ZXY)" << "\n";
	std::cout << "-f Transformfilename \t File with transformation parameter" << "\n";
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

	// input handling
	std::vector<double> Translation {0,0,0};
	std::vector<double> Rotation {0,0,0};
    std::string OutputFilename{ "" };
    std::string TransformationFileName{ "" };
	for (auto i = 1; i < argc; ++i)
	{
		if( std::string(argv[i]) == "-o")
		{
			++i;
			OutputFilename = argv[i];
		}
		if (std::string(argv[i]) == "-f")
		{
			++i;
			TransformationFileName = argv[i];
		}
		if (std::string(argv[i]) == "-t")
		{
			++i;
			Translation[0] = std::atof(argv[i]);
			++i;
			Translation[1] = std::atof(argv[i]);
			++i;
			Translation[2] = std::atof(argv[i]);
		}
		if (std::string(argv[i]) == "-r")
		{
			++i;
			Rotation[0] = std::atof(argv[i]);
			++i;
			Rotation[1] = std::atof(argv[i]);
			++i;
			Rotation[2] = std::atof(argv[i]);
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


    EulerTransformPointer T = nullptr;
    BSplineTransformPointer B = nullptr;
	if(!TransformationFileName.empty())
	{
		auto TransformFileReader = TransformReaderType::New();
		TransformFileReader->SetFileName(TransformationFileName);
		try
		{
			TransformFileReader->Update();
		}
		catch (itk::ExceptionObject &EO)
		{
			EO.Print(std::cout);
			return EXIT_FAILURE;
		}

		const TransformReaderType::TransformListType * transforms =	TransformFileReader->GetTransformList();
		auto it = transforms->begin();
		if (!strcmp((*it)->GetNameOfClass(), "Euler3DTransform"))
		{
			T = dynamic_cast< EulerTransformType* >((*it).GetPointer());
		}
        if (!strcmp((*it)->GetNameOfClass(), "BSplineTransform"))
        {
            B = dynamic_cast<BSplineTransformType*>((*it).GetPointer());
        }
		T.Print(std::cout);
	}
	else
	{
        T = EulerTransformType::New();
		T->SetRotation(deg2rad(Rotation[0]), deg2rad(Rotation[1]), deg2rad(Rotation[2]));
		EulerTransformType::OutputVectorType TranslationVector;
		for (size_t i = 0; i < Translation.size(); ++i)
		{
			TranslationVector[i] = Translation[i];
		}
		T->SetTranslation(TranslationVector);

		// Rotate over image center before translate
		auto origin = ImageReader->GetOutput()->GetOrigin();
		auto spacing = ImageReader->GetOutput()->GetSpacing();
		auto size = ImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
		origin[0] += spacing[0] * static_cast<double>(size[0]) / 2.0;
		origin[1] += spacing[1] * static_cast<double>(size[1]) / 2.0;
		origin[2] += spacing[2] * static_cast<double>(size[2]) / 2.0;
		T->SetCenter(origin);
	}

	auto TransformWriter = itk::TransformFileWriterTemplate<double>::New();
    (T != nullptr) ? TransformWriter->SetInput(T) : TransformWriter->SetInput(B); // writes only BSpline or Euler3D
	const std::experimental::filesystem::path Filename(OutputFilename);
	TransformWriter->SetFileName(Filename.stem().string() + ".tfm");
    try
    {
        std::cout << Filename.stem().string() + ".tfm" << std::endl;
        TransformWriter->Update();
    }
    catch (itk::ExceptionObject &E)
    {
        E.Print(std::cout);
    }
	

	auto ResampleFilter = ResampleImageFilterType::New();
	ResampleFilter->SetInput(ImageReader->GetOutput());
	ResampleFilter->SetDefaultPixelValue(0);

	ResampleFilter->SetSize(ImageReader->GetOutput()->GetLargestPossibleRegion().GetSize());
	ResampleFilter->SetOutputSpacing(ImageReader->GetOutput()->GetSpacing());

    if(T != nullptr)
    {
        ResampleFilter->SetTransform(T);
        const auto OutputOrigin = T->GetInverseTransform()->TransformPoint(ImageReader->GetOutput()->GetOrigin());
        ResampleFilter->SetOutputOrigin(OutputOrigin);
        const auto InverseTransform = T->Clone();
        T->GetInverse(InverseTransform);
        const auto OutputDirection = InverseTransform->GetMatrix() * ImageReader->GetOutput()->GetDirection();
        ResampleFilter->SetOutputDirection(OutputDirection);
    }

    if(B != nullptr)
    {
        ResampleFilter->SetTransform(B);
        ResampleFilter->SetOutputOrigin(ImageReader->GetOutput()->GetOrigin());
        ResampleFilter->SetOutputDirection(ImageReader->GetOutput()->GetDirection());
    }


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
	ImageWriter->SetFileName(OutputFilename);
	
	try
	{
        std::cout << OutputFilename << std::endl;
		ImageWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
	
}
