#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "rtkThreeDCircularProjectionGeometryXMLFileReader.h"
#include "itkRTKForwardProjectionImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkEuler3DTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"

using PixelType = float;
static constexpr unsigned int Dim2 = 2;
static constexpr unsigned int Dim3 = 3;
using Image2DType = itk::CudaImage<PixelType, Dim2>;
using Image3DType = itk::CudaImage<PixelType, Dim3>;
using Image2DReaderType = itk::ImageFileReader<Image2DType>;
using Image3DReaderType = itk::ImageFileReader<Image3DType>;
using GeometryReaderType = rtk::ThreeDCircularProjectionGeometryXMLFileReader;
using InterpolatorType = itk::LinearInterpolateImageFunction<Image2DType>;

using ProjectionFilterType = itk::RTKForwardProjectionImageFilter<Image3DType, Image3DType>;
using ExtractImagFilterType = itk::ExtractImageFilter<Image3DType, Image2DType>;
using Euler3DTransformType = itk::Euler3DTransform<double>;

using MSMetricType = itk::MeanSquaresImageToImageMetric<Image2DType, Image2DType>;
using IdentityTransformType = itk::IdentityTransform<double, Dim2>;

void
PrintUsage(std::string programme)
{
	std::cout << "USAGE:\n";
	std::cout << programme << " Fixed2DImageFile Moving3DImageFile ProjectionGeometry <options> \n";
	std::cout << "\nOptions:\n";
	std::cout << "-o OutputFilename" << "\t\t\t" << "name for output file\n";
	std::cout << "-x x_min x_max" << "\t\t\t" << "range of translation in x direction\n";
	std::cout << "-y y_min y_max" << "\t\t\t" << "range of translation in y direction\n";
	std::cout << "-z z_min z_max" << "\t\t\t" << "range of translation in z direction\n";
	std::cout << "-s x y z" << "\t\t\t" << "stepsize in x y z direction\n";

	std::cout << "\nInformation:\n";
	std::cout << "The default range is [-100,100] with a step size of 1. No default output.\n";
	std::cout << std::endl;
}

template<class ImageType>
void WriteImage(typename ImageType::Pointer image, Euler3DTransformType::ParametersType params, std::string name = "noname");


int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		PrintUsage(argv[0]);
		return EXIT_FAILURE;
	}

	// Default settings
	std::string OutputFilename("");
	std::pair<double, double> rangeX = { -100,100 };
	std::pair<double, double> rangeY = { -100,100 };
	std::pair<double, double> rangeZ = { -100,100 };
	std::vector<double> stepsize = { 1,1,1 };

	// Handling the input
	for (auto i = 4; i < argc; ++i)
	{
		if (std::string(argv[i]) == "-o")
		{
			++i;
			OutputFilename = argv[i];
		}
		if (std::string(argv[i]) == "-x")
		{
			++i;
			rangeX.first = std::stod(argv[i]);
			++i;
			rangeX.second = std::stod(argv[i]);
		}
		if (std::string(argv[i]) == "-y")
		{
			++i;
			rangeY.first = std::stod(argv[i]);
			++i;
			rangeY.second = std::stod(argv[i]);
		}
		if (std::string(argv[i]) == "-z")
		{
			++i;
			rangeZ.first = std::stod(argv[i]);
			++i;
			rangeZ.second = std::stod(argv[i]);
		}
		if (std::string(argv[i]) == "-s")
		{
			++i;
			stepsize[0] = std::stod(argv[i]);
			++i;
			stepsize[1] = std::stod(argv[i]);
			++i;
			stepsize[2] = std::stod(argv[i]);
		}

	}

	// Print state
	std::cout << "Fixed ImageFile: \t\t " << argv[1] << "\n";
	std::cout << "MovingImageFile: \t\t " << argv[2] << "\n";
	std::cout << "ProjectionGeometryFile: \t " << argv[3] << "\n";
	std::cout << "Range X: \t\t\t " << rangeX.first << " to " << rangeX.second << "\n";
	std::cout << "Range Y: \t\t\t " << rangeY.first << " to " << rangeY.second << "\n";
	std::cout << "Range Z: \t\t\t " << rangeZ.first << " to " << rangeZ.second << "\n";
	std::cout << "Stepsize [x,y,z]: \t\t [" << stepsize[0] << "," << stepsize[1] << "," << stepsize[2] << "]\n";
	std::cout << "OutputFile: \t\t\t " << OutputFilename << "\n";
	std::cout << std::endl;

	// Read fixed 2D image file
	auto FixedFileReader = Image2DReaderType::New();
	FixedFileReader->SetFileName(argv[1]);
	auto FixedImage = FixedFileReader->GetOutput();
	// Read moving 3D image file
	auto MovingFileReader = Image3DReaderType::New();
	MovingFileReader->SetFileName(argv[2]);
	auto MovingImage = MovingFileReader->GetOutput();

	// Read Projection geometry file
	auto GeometryReader = GeometryReaderType::New();
	GeometryReader->SetFilename(argv[3]);
	GeometryReader->GenerateOutputInformation();
	auto Geometry = GeometryReader->GetOutputObject();
	Geometry->Print(std::cout);

	try
	{
		FixedImage->Update();
		MovingImage->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	// Create projection components
	auto Transform = Euler3DTransformType::New();
	Transform->SetIdentity();

	auto ProjectionFilter = ProjectionFilterType::New();
	ProjectionFilter->SetInput(1, MovingImage);
	ProjectionFilter->SetDetectorImage<Image2DType>(FixedImage);
	ProjectionFilter->SetTransform(Transform.GetPointer());
	ProjectionFilter->SetBaseGeometry(GeometryReader->GetOutputObject());
	ProjectionFilter->InPlaceOff();

	auto ExtractImageFilter = ExtractImagFilterType::New();
	ExtractImageFilter->SetDirectionCollapseToSubmatrix();
	Image3DType::RegionType ExtractionRegion;
	Image3DType::IndexType ExtractionIndex;
	ExtractionIndex.Fill(0);
	ExtractionRegion.SetIndex(ExtractionIndex);
	Image3DType::SizeType ExtractionSize;
	ExtractionSize[0] = FixedImage->GetLargestPossibleRegion().GetSize()[0];
	ExtractionSize[1] = FixedImage->GetLargestPossibleRegion().GetSize()[1];
	ExtractionSize[2] = 0;
	ExtractionRegion.SetSize(ExtractionSize);
	ExtractImageFilter->SetExtractionRegion(ExtractionRegion);
	ExtractImageFilter->SetInput(ProjectionFilter->GetOutput());
	try
	{
		ExtractImageFilter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	//WriteImage<Image2DType>(ExtractImageFilter->GetOutput(), Transform->GetParameters(), "first");

	// Set up metric
	auto IdentityTransform = IdentityTransformType::New();
	IdentityTransform->SetIdentity();

	auto Interpolator = InterpolatorType::New();

	auto MSMetric = MSMetricType::New();
	MSMetric->SetInterpolator(Interpolator);
	MSMetric->SetTransform(IdentityTransform);
	MSMetric->SetFixedImage(FixedImage);
	MSMetric->SetFixedImageRegion(FixedImage->GetLargestPossibleRegion());
	MSMetric->SetMovingImage(ExtractImageFilter->GetOutput());

	Euler3DTransformType::ParametersType params(Transform->GetNumberOfParameters());
	params.Fill(0.0);

	try
	{
		MSMetric->Initialize();
	}
	catch(itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}
	

	std::ofstream OutputStream;
	auto ToFile = false;

	if (OutputFilename.compare("") != 0)
	{
		OutputStream.open(OutputFilename, std::ios::trunc);
		if (!OutputStream.is_open()) {
			std::cout << "File could not be opened." << std::endl;
		}
		ToFile = true;
	}

	for (auto x = rangeX.first; x <= rangeX.second; x += stepsize[0])
	{
		params(3) = x;
		for (auto y = rangeY.first; y <= rangeY.second; y += stepsize[1])
		{
			params(4) = y;
			for (auto z = rangeZ.first; z <= rangeZ.second; z += stepsize[2])
			{
				params(5) = z;

				try
				{
					Transform->SetParameters(params);
					ProjectionFilter->Modified();
					ProjectionFilter->Update();
					ExtractImageFilter->UpdateLargestPossibleRegion();
					MSMetric->GetMovingImage()->GetSource()->Update();
					auto value = MSMetric->GetValue(IdentityTransform->GetParameters());
					
					if (ToFile)
					{
						OutputStream << MSMetric.GetPointer() << " ";
						OutputStream << value << " ";
						OutputStream << Transform->GetParameters();
						OutputStream << std::endl;
					}
					else
					{
						std::cout << Transform->GetParameters() << std::endl;
					}

					//WriteImage<Image2DType>(ExtractImageFilter->GetOutput(), Transform->GetParameters());
				}
				catch (itk::ExceptionObject &EO)
				{
					EO.Print(std::cout);
					if (ToFile)
					{
						OutputStream << MSMetric.GetPointer() << " ";
						OutputStream << "NaN" << " ";
						OutputStream << Transform->GetParameters();
						OutputStream << std::endl;
					}
					else
					{
						std::cout << Transform->GetParameters() << std::endl;
					}
				}
			}
		}
	}

	OutputStream.close();

	return EXIT_SUCCESS;


}

template<class ImageType>
void
WriteImage(typename ImageType::Pointer image, Euler3DTransformType::ParametersType params, std::string name)
{
	auto Writer = itk::ImageFileWriter<ImageType>::New();
	Writer->SetInput(image);

	std::string OutputFilename("MI");
	if(name.compare("noname") != 0)
	{
		OutputFilename = name;
	}
	else
	{
		OutputFilename.append("_");
		OutputFilename.append(std::to_string(params(3)));
		OutputFilename.append("_");
		OutputFilename.append(std::to_string(params(4)));
		OutputFilename.append("_");
		OutputFilename.append(std::to_string(params(5)));
	}
	OutputFilename.append(".nrrd");
	Writer->SetFileName(OutputFilename);

	try
	{
		Writer->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
	}
}