#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "gdcmUIDGenerator.h"
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
#include "itkImageSeriesWriter.h"

#include "itkCommand.h"
#include "itkVersion.h"

// ITK typedef
using PixelType = float;
static constexpr int Dimension = 3;
using ImageType = itk::Image<PixelType, Dimension>;

using SpatialImageObjectType = itk::ImageSpatialObject<Dimension, PixelType>;

using OutputPixelType = unsigned int;
static constexpr int OutputDimension = 2;
using OutputImageType = itk::Image<OutputPixelType, OutputDimension>;

using ImageSeriesReaderType = itk::ImageSeriesReader<ImageType>;
using ImageSeriesWriterType = itk::ImageSeriesWriter<ImageType, OutputImageType>;
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
PrintUsage
(char* programname)
{
	std::cout << "USAGE: " << std::endl;
	std::cout << programname << " FixedDicomDirectory MovingDicomDirectory OutputDirectory" << std::endl;
}

ImageType::Pointer
Load3DVolume
(std::string, GDCMIOType::Pointer);

void
CopyDictionary
(itk::MetaDataDictionary&, itk::MetaDataDictionary&);

void
WriteTransformedVolume
(std::string, std::string, ImageType::Pointer, ImageSeriesReaderType::DictionaryArrayType, GDCMIOType::Pointer);

ImageSeriesReaderType::DictionaryArrayType
LoadDicomMetaDataDictionary
(const char*, ImageType::Pointer, GDCMIOType::Pointer);


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
	Resampler->SetDefaultPixelValue(-1024);

	try
	{
		Resampler->Update();
	}
	catch (itk::ExceptionObject EO)
	{
		std::cerr << "Resample error" << std::endl;
		EO.Print(std::cerr);
		return EXIT_FAILURE;
	}

	// Output dicom creation
	auto OutputDicomHeader = LoadDicomMetaDataDictionary(argv[1], Resampler->GetOutput(), FixedGDCMIO);
	WriteTransformedVolume(argv[3], argv[1], Resampler->GetOutput(), OutputDicomHeader, FixedGDCMIO);
	
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

void
CopyDictionary
(itk::MetaDataDictionary &SrcDict, itk::MetaDataDictionary &DestDict)
{
	typedef itk::MetaDataDictionary DictionaryType;

	DictionaryType::ConstIterator Iter = SrcDict.Begin();
	DictionaryType::ConstIterator End = SrcDict.End();
	typedef itk::MetaDataObject< std::string > MetaDataStringType;

	while (Iter != End)
	{
		auto Entry = Iter->second;

		MetaDataStringType::Pointer EntryValue = dynamic_cast<MetaDataStringType *>(Entry.GetPointer());
		if (EntryValue)
		{
			std::string TagKey = Iter->first;
			std::string TagValue = EntryValue->GetMetaDataObjectValue();
			itk::EncapsulateMetaData<std::string>(DestDict, TagKey, TagValue);
		}
		++Iter;
	}
}

void
WriteTransformedVolume
(std::string OutputDirectory, std::string InputDirectory, ImageType::Pointer Image, ImageSeriesReaderType::DictionaryArrayType Dict, GDCMIOType::Pointer GDCMIO)
{
	// setup the directories for the output
	itksys::SystemTools::MakeDirectory(OutputDirectory);

	auto FileNameGenerator = NamesGeneratorType::New();
	FileNameGenerator->SetInputDirectory(InputDirectory);
	auto FileNames = FileNameGenerator->GetInputFileNames();

	printf("Writing output series\n");
	auto SeriesWriter = ImageSeriesWriterType::New();
	SeriesWriter->SetInput(Image);
	SeriesWriter->SetImageIO(GDCMIO);
	FileNameGenerator->SetOutputDirectory(OutputDirectory);
	SeriesWriter->SetFileNames(FileNameGenerator->GetOutputFileNames());
	SeriesWriter->SetMetaDataDictionaryArray(&Dict);

	try
	{
		SeriesWriter->Update();
		std::cout << "\tDONE.\n" << std::endl;
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cerr << "Exception thrown while writing the series " << std::endl;
		std::cerr << EO << std::endl;
	}
}

ImageSeriesReaderType::DictionaryArrayType
LoadDicomMetaDataDictionary
(const char* Directory, ImageType::Pointer OutputImage, GDCMIOType::Pointer GDCMIO)
{
	std::cout << "Copy and update header data..." << std::endl;
	auto FileNameGenerator = NamesGeneratorType::New();

	// reading the volume meta data
	FileNameGenerator->SetInputDirectory(Directory);
	const ImageSeriesReaderType::FileNamesContainer & FileNames = FileNameGenerator->GetInputFileNames();

	auto Reader = ImageSeriesReaderType::New();
	Reader->SetImageIO(GDCMIO);
	Reader->SetFileNames(FileNames);
	Reader->Update();
	ImageSeriesReaderType::DictionaryArrayType OutputArray;

	auto InputDict = (*(Reader->GetMetaDataDictionaryArray()))[0];

	// Generate new SeriesID and Frame of reference (they all share the same spatial reference)
	gdcm::UIDGenerator SUID;
	std::string SeriesUID = SUID.Generate();
	gdcm::UIDGenerator FUID;
	std::string FrameOfReferenceUID = FUID.Generate();
	// Retireve the Study ID and the SOP Class ID
	std::string StudyUID;
	std::string SOPClassUID;
	itk::ExposeMetaData<std::string>(*InputDict, "0020|000d", StudyUID);
	itk::ExposeMetaData<std::string>(*InputDict, "0008|0016", SOPClassUID);
	GDCMIO->KeepOriginalUIDOn();

	for (unsigned int f = 0; f < Reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2]; f++)
	{
		// Create a new dictionary for this slice
		ImageSeriesReaderType::DictionaryRawPointer Dict = new ImageSeriesReaderType::DictionaryType;

		// Copy the dictionary from the first slice
		CopyDictionary(*InputDict, *Dict);

		// Set the UID's for the study, series, SOP  and frame of reference
		itk::EncapsulateMetaData<std::string>(*Dict, "0020|000d", StudyUID);
		itk::EncapsulateMetaData<std::string>(*Dict, "0020|000e", SeriesUID);
		itk::EncapsulateMetaData<std::string>(*Dict, "0020|0052", FrameOfReferenceUID);

		// Dicom share a SOPClass but every dicom has its own SOPInstanceUID and this is usually also the Media Storage SOPInstanceUID
		gdcm::UIDGenerator SOPUID;
		std::string SOPInstanceUID = SOPUID.Generate();
		itk::EncapsulateMetaData<std::string>(*Dict, "0008|0018", SOPInstanceUID);
		itk::EncapsulateMetaData<std::string>(*Dict, "0002|0003", SOPInstanceUID);

		// Change fields that are slice specific
		std::ostringstream Value;
		Value.str("");
		Value << f + 1;

		// Image Number
		itk::EncapsulateMetaData<std::string>(*Dict, "0020|0013", Value.str());

		// Series Description - Append new description to current series description
		std::string OldSeriesDesc;
		itk::ExposeMetaData<std::string>(*InputDict, "0008|103e", OldSeriesDesc);

		Value.str("");
		Value << OldSeriesDesc << ": + TRANSFORMED";

		// This is an long string and there is a 64 character limit in the standard
		unsigned LengthDesc = Value.str().length();

		std::string SeriesDesc(Value.str(), 0, LengthDesc > 64 ? 64 : LengthDesc);
		itk::EncapsulateMetaData<std::string>(*Dict, "0008|103e", SeriesDesc);

		// Series Number can be any number, Should be the same for all objects in study
		Value.str("");
		Value << 1001;
		itk::EncapsulateMetaData<std::string>(*Dict, "0020|0011", Value.str());

		// Derivation description - how this image was derived
		Value.str("");
		Value << "ITK: " << ITK_SOURCE_VERSION;

		LengthDesc = Value.str().length();
		std::string DerivationDesc(Value.str(), 0, LengthDesc > 1024 ? 1024 : LengthDesc);
		itk::EncapsulateMetaData<std::string>(*Dict, "0008|2111", DerivationDesc);

		// Image Position Patient: This is calculated by computing the
		// physical coordinate of the first pixel in each slice.
		ImageType::PointType Position;
		ImageType::IndexType Index;
		Index[0] = 0;
		Index[1] = 0;
		Index[2] = f;
		OutputImage->TransformIndexToPhysicalPoint(Index, Position);

		Value.str("");
		Value << Position[0] << "\\" << Position[1] << "\\" << Position[2];
		itk::EncapsulateMetaData<std::string>(*Dict, "0020|0032", Value.str());
		// Slice Location: For now, we store the z component of the Image Position Patient.
		Value.str("");
		Value << Position[2];
		itk::EncapsulateMetaData<std::string>(*Dict, "0020|1041", Value.str());

		// Save the dictionary
		OutputArray.push_back(Dict);
	}
	std::cout << "\tDONE.\n" << std::endl;
	return OutputArray;
}