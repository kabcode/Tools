#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"

#include "itkImageRegistrationMethodv4.h"
#include "itkCenteredTransformInitializer.h"
#include "itkEuler3DTransform.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkGradientDescentOptimizerv4.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageSpatialObject.h"

#include "itkResampleImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkCommand.h"

// ITK typedef
using PixelType = float;
static constexpr int Dimension = 3;
using ImageType = itk::Image<PixelType, Dimension>;

using SpatialImageObjectType = itk::ImageSpatialObject<Dimension, PixelType>;

using OutputPixelType = unsigned int;
static constexpr int OutputDimension = 2;
using OutputImageType = itk::Image<OutputPixelType, OutputDimension>;

using ImageSeriesReaderType = itk::ImageSeriesReader<ImageType>;
using ImageWriterType = itk::ImageFileWriter<ImageType>;

// GDCM typedefs
using GDCMIOType = itk::GDCMImageIO;
using NamesGeneratorType = itk::GDCMSeriesFileNames;

// Registration typedefs
using EulerTransformType = itk::Euler3DTransform<double>;
using CenteredTransfromInitializerType = itk::CenteredTransformInitializer<EulerTransformType, ImageType, ImageType>;
using MIMetricType = itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType>;
using ScalesEstimatorType = itk::RegistrationParameterScalesFromPhysicalShift<MIMetricType>;
using GDOptimizerType = itk::GradientDescentOptimizerv4;
using LinearInterpolatorType = itk::LinearInterpolateImageFunction<ImageType>;
using RegistrationType = itk::ImageRegistrationMethodv4<ImageType, ImageType, EulerTransformType>;

using ThresholdImageFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;

// Resampling typedefs
using ResamplerType = itk::ResampleImageFilter<ImageType, ImageType>;


// Command/Observer
class CommandIterationUpdate : public itk::Command
{
public:
	using Self = CommandIterationUpdate;
	using Superclass = itk::Command;
	using Pointer = itk::SmartPointer<Self>;
	itkNewMacro(Self);
protected:
	CommandIterationUpdate() = default;
public:
	using OptimizerPointer = const GDOptimizerType *;
	void Execute(itk::Object *caller, const itk::EventObject & event) override
	{
		Execute((const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event) override
	{
		auto optimizer = static_cast< OptimizerPointer >(object);
		if (!itk::IterationEvent().CheckEvent(&event))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << std::endl;
	}
};

// Helper functions
void
PrintUsage(char* programname)
{
	std::cout << "USAGE: " << std::endl;
	std::cout << programname << " FixedDicomDirectory MovingDicomDirectory OutputDirectory" << std::endl;
}

ImageType::Pointer
Load3DVolume
(std::string, GDCMIOType::Pointer);

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		PrintUsage(argv[0]);
		return EXIT_FAILURE;
	}


	// Read fixed DICOM series
	auto FixedGDCMIO = GDCMIOType::New();
	auto FixedImage = Load3DVolume(argv[1], FixedGDCMIO);

	// Read moving DICOM series
	auto MovingGDCMIO = GDCMIOType::New();
	auto MovingImage = Load3DVolume(argv[2], MovingGDCMIO);

	// Create fixed image mask
	auto ThresholdImageFilter = ThresholdImageFilterType::New();
	ThresholdImageFilter->SetInput(FixedImage);
	ThresholdImageFilter->SetUpperThreshold(300);
	ThresholdImageFilter->SetLowerThreshold(-600);
	auto MaskedImage = ThresholdImageFilter->GetOutput();

	auto MaskWriter = ImageWriterType::New();
	MaskWriter->SetFileName("Mask.nrrd");
	MaskWriter->SetInput(MaskedImage);
	
	try
	{
		MaskedImage->Update();
		MaskWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cerr << "Thresholding error" << std::endl;
		EO.Print(std::cerr);
		return EXIT_FAILURE;
	}

	auto FixedImageMask = SpatialImageObjectType::New();
	FixedImageMask->SetImage(MaskedImage);
	
	// Define Registration components
	auto RegistrationMethod = RegistrationType::New();
	auto EulerTransform = EulerTransformType::New();
	auto MutualInformationMetric = MIMetricType::New();
	auto ScalesEstimator = ScalesEstimatorType::New();
	auto GradientDescentOptimizer = GDOptimizerType::New();
	auto LinearInterpolator = LinearInterpolatorType::New();

	// Initial transform of the images
	auto InitialTransform = CenteredTransfromInitializerType::New();
	InitialTransform->SetFixedImage(FixedImage);
	InitialTransform->SetMovingImage(MovingImage);
	InitialTransform->SetTransform(EulerTransform);
	InitialTransform->GeometryOn();
	InitialTransform->InitializeTransform();

	// Setup registration components
	MutualInformationMetric->SetNumberOfHistogramBins(50);
	MutualInformationMetric->SetFixedImageMask(FixedImageMask);
	ScalesEstimator->SetMetric(MutualInformationMetric);
	ScalesEstimator->SetTransformForward(true);
	ScalesEstimator->SetSmallParameterVariation(1.0);
	GradientDescentOptimizer->SetLearningRate(1.0);
	GradientDescentOptimizer->SetNumberOfIterations(200);
	GradientDescentOptimizer->SetMinimumConvergenceValue(1e-6);
	GradientDescentOptimizer->SetConvergenceWindowSize(10);
	GradientDescentOptimizer->SetScalesEstimator(ScalesEstimator);

	// Link components to registration
	RegistrationMethod->SetMetric(MutualInformationMetric);
	RegistrationMethod->SetMetricSamplingStrategy(RegistrationType::MetricSamplingStrategyType::RANDOM);
	RegistrationMethod->SetMetricSamplingPercentage(0.01);
	RegistrationMethod->SetInitialTransform(EulerTransform);
	RegistrationMethod->SetOptimizer(GradientDescentOptimizer);
	
	// Set level resolution parameters
	int numberoflevels = 3;
	RegistrationMethod->SetNumberOfLevels(numberoflevels);
	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(numberoflevels);
	shrinkFactorsPerLevel[0] = 4;
	shrinkFactorsPerLevel[1] = 2;
	shrinkFactorsPerLevel[2] = 1;
	RegistrationMethod->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);

	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(numberoflevels);
	smoothingSigmasPerLevel[0] = 2;
	smoothingSigmasPerLevel[1] = 1;
	smoothingSigmasPerLevel[2] = 0;
	RegistrationMethod->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	RegistrationMethod->SmoothingSigmasAreSpecifiedInPhysicalUnitsOn();

	// Set images
	RegistrationMethod->SetFixedImage(FixedImage);
	RegistrationMethod->SetMovingImage(MovingImage);

	// Add observer
	auto OptimizerObserver = CommandIterationUpdate::New();
	GradientDescentOptimizer->AddObserver(itk::IterationEvent(), OptimizerObserver);

	// Run registration
	try
	{
		RegistrationMethod->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cerr << "Registration error" << std::endl;
		EO.Print(std::cerr);
		return EXIT_FAILURE;
	}

	// Resample moving image to output image
	auto Resampler = ResamplerType::New();
	Resampler->SetInput(MovingImage);
	Resampler->SetReferenceImage(FixedImage);
	Resampler->UseReferenceImageOn();
	Resampler->SetTransform(EulerTransform);

	auto ImageWriter = ImageWriterType::New();
	ImageWriter->SetFileName("Out.nrrd");
	ImageWriter->SetInput(Resampler->GetOutput());

	try
	{
		ImageWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cerr << "Writing error" << std::endl;
		EO.Print(std::cerr);
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;

}

ImageType::Pointer
Load3DVolume
(std::string Directory, GDCMIOType::Pointer GDCMIO)
{
	std::cout << "Load Volume from: " << Directory << std::endl;
	auto FileNameGenerator = NamesGeneratorType::New();

	// Reading the volume data
	FileNameGenerator->SetInputDirectory(Directory);
	auto FileNames = FileNameGenerator->GetInputFileNames();
	unsigned int NumberOfFilenames = FileNames.size();

	auto Reader = ImageSeriesReaderType::New();
	Reader->SetImageIO(GDCMIO);
	Reader->SetFileNames(FileNames);
	try
	{
		Reader->Update();
		std::cout << "\tNumber of files: " << NumberOfFilenames << "\n";
		std::cout << "\tDONE.\n" << std::endl;
		return Reader->GetOutput();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cerr << "Exception thrown while reading the image" << std::endl;
		std::cerr << EO << std::endl;
		return nullptr;
	}
}