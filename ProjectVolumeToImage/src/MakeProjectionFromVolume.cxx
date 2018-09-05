// STL includes
#include <filesystem>

// ITK includes
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"
#include "itkEuler3DTransform.h"
#include "itkResampleImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "gdcmUUIDGenerator.h"

// MORF includes
#include "itkRTKForwardProjectionImageFilter.h"
#include "rtkThreeDCircularProjectionGeometryXMLFile.h"

namespace FS = std::experimental::filesystem;
using UUIDGeneratorType = gdcm::UUIDGenerator;

using PixelType = float;
const unsigned int Dim3D = 3;
const unsigned int Dim2D = 2;
using InputImageType = itk::CudaImage<PixelType, Dim3D>;
using InternalImageType = itk::CudaImage<PixelType, Dim3D>;
using OutputImageType = itk::CudaImage<PixelType, Dim2D>;
using InternalImageIterator = itk::ImageRegionIterator<InternalImageType>;
using InternalImageConstIterator = itk::ImageRegionConstIterator<InternalImageType>;

using ImageReaderType = itk::ImageFileReader<InputImageType>;
using ImageWriterType = itk::ImageFileWriter<InternalImageType>;
using TransformWriterType = itk::TransformFileWriter;
using OutputImageWriterType = itk::ImageFileWriter<OutputImageType>;
using EulerTransformType = itk::Euler3DTransform<double>;
using ResampleImageFilterType = itk::ResampleImageFilter<InternalImageType, InternalImageType>;
using ExtractionFilterType = itk::ExtractImageFilter<InternalImageType, OutputImageType>;
using ImageInformationChangeFilterType = itk::ChangeInformationImageFilter<InternalImageType>;

using ProjectorType = itk::RTKForwardProjectionImageFilter<InputImageType, InternalImageType>;
using ProjectionGeometryType = rtk::ThreeDCircularProjectionGeometry;
using ProjectionGeometryWriterType = rtk::ThreeDCircularProjectionGeometryXMLFileWriter;

const double Pi = itk::Math::pi;

void
PrintUsage(char * programmname)
{
	std::cout << "\nInformation to " << programmname << std::endl;
	std::cout << "Usage: " << std::endl;
	std::cout << programmname << " <options> <inputoptions>" << std::endl;
	std::cout << std::endl;
	std::cout << "The default position of the geometry is the source pointing towards the positive y-axis (at patient in HFS it's along anterior-posterior).\n";
	std::cout << "The default values are SSD 1200mm, SID 800mm, spacing 0.15mm, size 1200x1400 px, isocenter at (0,0,0).\n";
	std::cout << "All values are in mm or degree units and carried out in patient space." << std::endl;
	std::cout << std::endl;
	std::cout << "options:" << std::endl;
	std::cout << "-dox" << "\t\t\t" << "Detector offset towards x" << std::endl;
	std::cout << "-doz" << "\t\t\t" << "Detector offset towards z" << std::endl;
	std::cout << "-dsx" << "\t\t\t" << "Detector size in x direction" << std::endl;
	std::cout << "-dsz" << "\t\t\t" << "Detector size in z direction" << std::endl;
	std::cout << "-spx" << "\t\t\t" << "Detector spacing in x direction" << std::endl;
	std::cout << "-spz" << "\t\t\t" << "Detector spacing in z direction" << std::endl;
	std::cout << "-sdd" << "\t\t\t" << "Source Detector distance" << std::endl;
	std::cout << "-sid" << "\t\t\t" << "Source Isocenter distance" << std::endl;
	std::cout << "-sox" << "\t\t\t" << "Source offset towards x" << std::endl;
	std::cout << "-soz" << "\t\t\t" << "Source offset towards z" << std::endl;
	std::cout << "-px" << "\t\t\t" << "Translation of geometry towards x" << std::endl;
	std::cout << "-py" << "\t\t\t" << "Translation of geometry towards y" << std::endl;
	std::cout << "-pz" << "\t\t\t" << "Translation of geometry towards z" << std::endl;
	std::cout << "-ga" << "\t\t\t" << "Rotation around y axis (if patient HFS, Gantry angle)" << std::endl;
	std::cout << "-oa" << "\t\t\t" << "Rotation around x axis (patient axis in HFS, OutOfPlane angle)" << std::endl;
	std::cout << "-ia" << "\t\t\t" << "Rotation around z axis (patient axis in HFS, InPlane angle)" << std::endl;
	std::cout << "-g" << "\t\t\t" << "Provided geometry file" << std::endl;
	std::cout << "-src" << "\t\t\t" << "Write file containing a single pixel at source position" << std::endl;
	std::cout << "-h" << "\t\t\t" << "Print this information" << std::endl;

	std::cout << "\ninputoptions:" << std::endl;
	std::cout << "-d" << "\t\t\t" << "Input directory for DICOM" << std::endl;
	std::cout << "-f" << "\t\t\t" << "Input file for other formats" << std::endl;

	std::cout << "\nExamples:" << std::endl;
	std::cout << programmname << " -sdd 1200 -sox 10 -ga 45 -f Test.nrrd" <<std::endl;
	std::cout << programmname << " -g GeometryFile.xml -d path/to/dicom/directory" << std::endl;
}

inline float DegreesToRadians(const float degrees) { return degrees * Pi / 180; };

int TestProjection(std::string);

int main(int argc, char* argv[])
{
	// check input arguments
	if (argc < 3)
	{
		PrintUsage(argv[0]);
		return EXIT_FAILURE;
	}

	// control flags
	auto GeometryProvided = false;
	auto GeometryFile("");
	auto DICOMDirectory = false;
	auto ImageSource("");
	auto WriteSrcPxl = false;

	// Detector image resembles the actual detector
	auto DetectorImage  = InternalImageType::New();
	auto ProjectorImage = InternalImageType::New();
	InternalImageType::IndexType Index;
	Index.Fill(0);
	
	InternalImageType::SizeType Size;
	Size[0] = 1200;
	Size[1] = 1400;
	Size[2] = 1;
	
	InternalImageType::RegionType Region(Index, Size);
	DetectorImage->SetRegions(Region);
	DetectorImage->Allocate();

	InternalImageType::RegionType Region2(Index, Size);
	ProjectorImage->SetRegions(Region2);
	ProjectorImage->Allocate();

	InternalImageType::SpacingType Spacing;
	Spacing[0] = 0.15;
	Spacing[1] = 0.15;
	Spacing[2] = 1;
	DetectorImage->SetSpacing(Spacing);
	ProjectorImage->SetSpacing(Spacing);

	InternalImageType::PointType Origin;
	Origin.Fill(0);
	DetectorImage->SetOrigin(Origin);
	ProjectorImage->SetOrigin(Origin);

	// set default projection parameter in device space (z pointing from detector towards source)
	auto ProjectionGeometry = ProjectionGeometryType::New();
	auto SID = 800.;
	auto SDD = 1200.;

	InternalImageType::PointType SourcePosition;
	SourcePosition[0] = 0;
	SourcePosition[1] = 0;
	SourcePosition[2] = SID;

	InternalImageType::PointType DetectorPosition;
	DetectorPosition[0] = - static_cast<double>(Size[0]) / 2. * Spacing[0];
	DetectorPosition[1] = - static_cast<double>(Size[1]) / 2. * Spacing[1];
	DetectorPosition[2] = - (SDD - SID);

	ProjectionGeometryType::VectorType RowDirection;
	RowDirection[0] = 1;
	RowDirection[1] = 0;
	RowDirection[2] = 0;

	ProjectionGeometryType::VectorType ColumnDirection;
	ColumnDirection[0] = 0;
	ColumnDirection[1] = 1;
	ColumnDirection[2] = 0;

	// Initial transformation in patient space
	auto GantryAngle     = 0.; // patient y-axis
	auto InPlaneAngle    = 0.; // patient z-axis
	auto OutOfPlaneAngle = 0.; // patient x-axis
	auto ShiftX          = 0.;
	auto ShiftY          = 0.;
	auto ShiftZ          = 0.;

	// read user input parameter
	auto DetectorChanged = false;
	for (auto i = 1; i < argc; ++i)
	{
		// Detector base properties (size, spacing, sid, sdd)
		if (std::string(argv[i]).substr(0, 4) == "-spx")	// spacing in x
		{
			++i;
			Spacing[0] = atof(argv[i]);
			DetectorChanged = true;
		}
		if (std::string(argv[i]).substr(0, 4) == "-spy")	// spacing in z
		{
			++i;
			Spacing[1] = atof(argv[i]);
			DetectorChanged = true;
		}
		if (std::string(argv[i]).substr(0, 4) == "-dsx")	// size in x
		{
			++i;
			Size[0] = atof(argv[i]);
			DetectorChanged = true;
		}
		if (std::string(argv[i]).substr(0, 4) == "-dsy")	// size in z
		{
			++i;
			Size[1] = atof(argv[i]);
			DetectorChanged = true;
		}
		
		// Geometry advanced properties (offset detector, offset source)
		if (std::string(argv[i]).substr(0, 4) == "-sid")	// source isocenter distance
		{
			++i;
			SID = atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 4) == "-sdd")	// source detector distance
		{
			++i;
			SDD = atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 4) == "-dox")		// offset detector in x
		{
			++i;
			DetectorPosition[0] += atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 4) == "-doy")		// offset detector in z
		{
			++i;
			DetectorPosition[1] += atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 4) == "-sox")		// offset source in x
		{
			++i;
			SourcePosition[0] = atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 4) == "-soy")		// offset source in z
		{
			++i;
			SourcePosition[1] = atof(argv[i]);
		}

		// Projector movement properties
		if (std::string(argv[i]).substr(0, 3) == "-ga")		// gantry angle
		{
			++i;
			GantryAngle = atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 3) == "-ia")		// in-plane angle
		{
			++i;
			InPlaneAngle = atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 3) == "-oa")		// out-of-plane angle
		{
			++i;
			OutOfPlaneAngle = atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 3) == "-px")		// projector shift in x
		{
			++i;
			ShiftX = atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 3) == "-py")		// projector shift in y
		{
			++i;
			ShiftY = atof(argv[i]);
		}
		if (std::string(argv[i]).substr(0, 3) == "-pz")		// projector shift in z
		{
			++i;
			ShiftZ = atof(argv[i]);
		}

		if (std::string(argv[i]).substr(0, 2) == "-g")
		{
			GeometryProvided = true;
			++i;
			GeometryFile = argv[i];
		}
		if (std::string(argv[i]).substr(0, 4) == "-src")
		{
			WriteSrcPxl = true;
		}

		// general options
		if (std::string(argv[i]).substr(0,2 ) == "-h" )
		{
			PrintUsage(argv[0]);
			return EXIT_FAILURE;
		}

		// input options
		if (std::string(argv[i]).substr(0, 2) == "-f")
		{
			DICOMDirectory = false;
			++i;
			ImageSource = argv[i];
		}
		if (std::string(argv[i]).substr(0, 2) == "-d")
		{
			DICOMDirectory = true;
			++i;
			ImageSource = argv[i];
		}
	}

	// Adjust the default projector geometry to the user input
	SourcePosition[2]   = SID;
	DetectorPosition[2] = - (SDD - SID); // if SID or SDD changed

	if(DetectorChanged) // if Size or spacing has changed
	{
		Region.SetSize(Size);
		DetectorImage->SetRegions(Region);
		DetectorImage->Allocate();
		DetectorPosition[0] = -static_cast<double>(Size[0]) / 2. * Spacing[0];
		DetectorPosition[1] = -static_cast<double>(Size[1]) / 2. * Spacing[1];
		DetectorImage->SetSpacing(Spacing);

		Region2.SetSize(Size);
		ProjectorImage->SetRegions(Region2);
		ProjectorImage->Allocate();
		ProjectorImage->SetSpacing(Spacing);
	}

	// 1. Rotation R of projector geometry
	auto DetectorRotationTransform = EulerTransformType::New();
	DetectorRotationTransform->SetComputeZYX(true);
	DetectorRotationTransform->SetRotation(DegreesToRadians(OutOfPlaneAngle),
										   DegreesToRadians(InPlaneAngle),
									       DegreesToRadians(GantryAngle));
	const auto R = DetectorRotationTransform->GetMatrix();

	// 2. Translation of projector geometry
	auto ProjectorTranslationTransform = EulerTransformType::New();
	EulerTransformType::OutputVectorType Translation;
	Translation[0] = ShiftX;
	Translation[1] = ShiftY;
	Translation[2] = ShiftZ;
	ProjectorTranslationTransform->SetTranslation(Translation);

	// 3. Apply transformation matrices
	SourcePosition    = R * SourcePosition;
	DetectorPosition  = R * DetectorPosition;
	RowDirection      = R * RowDirection;
	ColumnDirection   = R * ColumnDirection;

	SourcePosition   = ProjectorTranslationTransform->TransformPoint(SourcePosition);
	DetectorPosition = ProjectorTranslationTransform->TransformPoint(DetectorPosition);

	// Add projection to geometry object
	if(GeometryProvided)
	{
		auto GeometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
		GeometryReader->SetFilename(GeometryFile);
		TRY_AND_EXIT_ON_ITK_EXCEPTION(GeometryReader->GenerateOutputInformation());

		ProjectionGeometry = GeometryReader->GetOutputObject();
	}
	else
	{
		ProjectionGeometry->AddProjection(SourcePosition, DetectorPosition, RowDirection, ColumnDirection);
	}
	
	// Read input image
	InputImageType::Pointer InputVolume;
	if(DICOMDirectory)
	{
		auto DCMIO = itk::GDCMImageIO::New();
		auto SeriesReader = itk::ImageSeriesReader<InputImageType>::New();
		SeriesReader->SetImageIO(DCMIO);
		auto NameGenerator = itk::GDCMSeriesFileNames::New();
		NameGenerator->SetDirectory(ImageSource);
		auto SeriesUID = NameGenerator->GetSeriesUIDs();
		std::string seriesIdentifier = SeriesUID.begin()->c_str();
		auto FileNames = NameGenerator->GetFileNames(seriesIdentifier);
		SeriesReader->SetFileNames(FileNames);
		try
		{
			SeriesReader->Update();
		}
		catch (...)
		{
			return EXIT_FAILURE;
		}
		InputVolume = SeriesReader->GetOutput();
	}
	else
	{
		auto ImageReader = itk::ImageFileReader<InputImageType>::New();
		ImageReader->SetFileName(ImageSource);
		try
		{
			ImageReader->Update();
		}
		catch (...)
		{
			return EXIT_FAILURE;
		}
		InputVolume = ImageReader->GetOutput();
	}
	
	try
	{
		InputVolume->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cout << "Load image error." << std::endl;
		EO.Print(std::cout);
	}

	// Set up projector
	auto Projector = ProjectorType::New();
	Projector->SetInput(0, ProjectorImage);
	Projector->SetInput(1, InputVolume);
	Projector->SetGeometry(ProjectionGeometry);

	try
	{
		Projector->Update();
	}
	catch(itk::ExceptionObject &EO)
	{
		std::cout << "Projection error" << std::endl;
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	// Copy ProjectionImage into DetectorImage buffer
	InternalImageConstIterator ProjectorImageIt(Projector->GetOutput(), Projector->GetOutput()->GetLargestPossibleRegion());
	ProjectorImageIt.GoToBegin();
	InternalImageIterator DetectorImageIt(DetectorImage, DetectorImage->GetLargestPossibleRegion());
	DetectorImageIt.GoToBegin();

	for (; !DetectorImageIt.IsAtEnd(); ++DetectorImageIt, ++ProjectorImageIt)
	{
		DetectorImageIt.Set(ProjectorImageIt.Get());
	}
	
	auto ChangeInformationFilter = ImageInformationChangeFilterType::New();
	ChangeInformationFilter->SetInput(DetectorImage);
	ChangeInformationFilter->SetOutputOrigin(DetectorPosition);
	ChangeInformationFilter->ChangeOriginOn();
	
	InternalImageType::DirectionType DetectorDirection;
	auto DetectorNormale = CrossProduct(RowDirection, ColumnDirection);
	for (auto i = 0; i < 3; ++i)
	{
		DetectorDirection[i][0] = RowDirection[i];
		DetectorDirection[i][1] = ColumnDirection[i];
		DetectorDirection[i][2] = DetectorNormale[i];
	}
	ChangeInformationFilter->SetOutputDirection(DetectorDirection);
	ChangeInformationFilter->ChangeDirectionOn();
	try
	{
		ChangeInformationFilter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cout << "Change image error" << std::endl;
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	// Check if subdirectory for projections exist otherwise create it
	auto InputDirectory = FS::current_path();
	FS::path ProjectionDirectory = InputDirectory.append("Projections");	
	if(!FS::exists(ProjectionDirectory))
	{
		if(!FS::create_directory(ProjectionDirectory))
		{
			std::cout << "Create projection directory error" << std::endl;
			return EXIT_FAILURE;
		}
	}

	// Write projection image and geometry
	auto ProjectionFileName = ProjectionDirectory;
	ProjectionFileName.append("P");
	auto DRRFileName = ProjectionDirectory;
	DRRFileName.append("D");
	auto GeometryFileName = ProjectionDirectory;
	GeometryFileName.append("G");
	UUIDGeneratorType UuidGenerator;
	std::string UUID = UuidGenerator.Generate();
	auto UUIDShort = UUID.substr(0, 7);
	ProjectionFileName += UUIDShort;
	ProjectionFileName += ".nrrd";
	DRRFileName        += UUIDShort;
	DRRFileName        += ".nrrd";
	GeometryFileName   += UUIDShort;
	GeometryFileName   +=".geom";
	std::cout << UUIDShort << std::endl;

	// src pixel
	if (WriteSrcPxl)
	{
		InternalImageType::IndexType SrcIdx;
		SrcIdx.Fill(0);
		InternalImageType::SizeType SrcSz;
		SrcSz.Fill(1);
		InternalImageType::RegionType SrcReg(SrcIdx, SrcSz);
		auto SrcPxl = InternalImageType::New();
		SrcPxl->SetRegions(SrcReg);
		SrcPxl->Allocate();
		SrcPxl->FillBuffer(1);
		SrcPxl->SetOrigin(SourcePosition);

		auto SrcPxlWriter = ImageWriterType::New();
		SrcPxlWriter->SetInput(SrcPxl);
		auto SrcPxlFileName = ProjectionDirectory;
		SrcPxlFileName.append("S");
		SrcPxlFileName += UUIDShort;
		SrcPxlFileName += ".nrrd";
		SrcPxlWriter->SetFileName(SrcPxlFileName.string());

		try
		{
			SrcPxlWriter->Update();
		}
		catch (itk::ExceptionObject &EO)
		{
			std::cout << "Write source point error." << std::endl;
			EO.Print(std::cout);
			return EXIT_FAILURE;
		}
	}

	auto ProjectionImageWriter = ImageWriterType::New();
	ProjectionImageWriter->SetInput(ChangeInformationFilter->GetOutput());
	ProjectionImageWriter->SetFileName(ProjectionFileName.string());

	try
	{
		ProjectionImageWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cout << "Write projection error." << std::endl;
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	auto ExtractionFilter = ExtractionFilterType::New();
	ExtractionFilter->SetInput(Projector->GetOutput());
	ExtractionFilter->SetDirectionCollapseToIdentity();
	auto inputRegion = Projector->GetOutput()->GetLargestPossibleRegion();
	auto size = inputRegion.GetSize();
	size[2] = 0; // we extract along y direction
	auto start = inputRegion.GetIndex();
	start[2] = 0;
	InternalImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);
	ExtractionFilter->SetExtractionRegion(desiredRegion);

	auto DRRImageWriter = OutputImageWriterType::New();
	DRRImageWriter->SetInput(ExtractionFilter->GetOutput());
	DRRImageWriter->SetFileName(DRRFileName.string());

	try
	{
		DRRImageWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cout << "DRR image write error." << std::endl;
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	// Write projection geometry to xml file
	auto ProjectionGeometryWriter = ProjectionGeometryWriterType::New();
	ProjectionGeometryWriter->SetObject(ProjectionGeometry);
	ProjectionGeometryWriter->SetFilename(GeometryFileName.string());
		
	try
	{
		ProjectionGeometryWriter->WriteFile();
	}
	catch (itk::ExceptionObject &EO)
	{
		std::cout << "Write geometry error." << std::endl;
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
