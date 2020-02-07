#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkJosephForwardProjectionImageFilter.h"
#include "itkRayCastInterpolateImageFunction.h"
#include "rtkZengForwardProjectionImageFilter.h"
#include "rtkCudaForwardProjectionImageFilter.h"

#include "rtkConstantImageSource.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "HighPrecisionTimer.h"

using OutputPixelType = float;
constexpr unsigned int Dimension = 3;
using OutputImageType = itk::CudaImage< OutputPixelType, Dimension >;

void WriteImageFile(OutputImageType::Pointer image, std::string& outputname, unsigned int size);

int main(int argc, char * argv[])
{

    rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
    geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
    geometryReader->SetFilename(argv[2]);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation())

        // Create a stack of empty projection images
        using ConstantImageSourceType = rtk::ConstantImageSource< OutputImageType >;
    ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
    constantImageSource->SetSize(OutputImageType::SizeType{ 512,512,1 });
    OutputImageType::SpacingType spacing;
    spacing.Fill(1);
    constantImageSource->SetSpacing(spacing);
    OutputImageType::PointType origin;
    origin[0] = -256;
    origin[1] = -256;
    origin[2] = 0;
    constantImageSource->SetOrigin(origin);

    // Adjust size according to geometry
    ConstantImageSourceType::SizeType sizeOutput;
    sizeOutput[0] = constantImageSource->GetSize()[0];
    sizeOutput[1] = constantImageSource->GetSize()[1];
    sizeOutput[2] = geometryReader->GetOutputObject()->GetGantryAngles().size();
    constantImageSource->SetSize(sizeOutput);

    // Input reader
    using ReaderType = itk::ImageFileReader<  OutputImageType >;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[1]);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->Update());

    auto ITKforwardProjection = itk::RayCastInterpolateImageFunction<OutputImageType>::New();
    auto focalpointRTK = geometryReader->GetOutputObject()->GetSourcePosition(0);
    itk::RayCastInterpolateImageFunction<OutputImageType>::InputPointType focalpoint;
    focalpoint[0] = focalpointRTK[0];
    focalpoint[1] = focalpointRTK[1] + 130;
    focalpoint[2] = focalpointRTK[2];

    std::cout << "Focal Point: " << focalpoint[0] << ", " << focalpoint[1] << ", " << focalpoint[2] << std::endl;

    ITKforwardProjection->SetFocalPoint(focalpoint);
    auto transform = itk::Euler3DTransform<>::New();
    transform->SetIdentity();
    ITKforwardProjection->SetTransform(transform);

    auto resampleImageFilter = itk::ResampleImageFilter<OutputImageType, OutputImageType>::New();
    resampleImageFilter->SetTransform(transform);
    resampleImageFilter->SetInput(reader->GetOutput());
    resampleImageFilter->SetInterpolator(ITKforwardProjection);
    resampleImageFilter->SetSize(sizeOutput);
    resampleImageFilter->SetOutputSpacing(spacing);
    auto ITKorigin = origin;
    ITKorigin[2] = -750;
    resampleImageFilter->SetOutputOrigin(ITKorigin);

    rtk::ForwardProjectionImageFilter<OutputImageType, OutputImageType>::Pointer forwardProjection;
    auto image = OutputImageType::New();


    for (unsigned int i = 1; i < 5; ++i)
    {
        std::string outputfilename{ "" };
        switch (i)
        {
        case(1):
            forwardProjection = rtk::JosephForwardProjectionImageFilter<OutputImageType, OutputImageType>::New();
            outputfilename.append("Joseph");
            break;
        case(2):
            forwardProjection = rtk::ZengForwardProjectionImageFilter<OutputImageType, OutputImageType>::New();
            outputfilename.append("Zeng");
            break;
        case(3):
            forwardProjection = rtk::CudaForwardProjectionImageFilter<OutputImageType, OutputImageType>::New();
            dynamic_cast<rtk::CudaForwardProjectionImageFilter<OutputImageType, OutputImageType>*>(forwardProjection.GetPointer())->SetStepSize(1);
            outputfilename.append("Cuda");
            break;
        case(4):

            outputfilename.append("ITK");
            break;
        default:
            std::cerr << "Unhandled --method value." << std::endl;
            return EXIT_FAILURE;
        }

        if (i == 1 || i == 2 || i == 3)
        {
            forwardProjection->SetInput(constantImageSource->GetOutput());
            forwardProjection->SetInput(1, reader->GetOutput());
            forwardProjection->SetGeometry(geometryReader->GetOutputObject());
            {
                auto Timer{ HighPrecisionTimer<TimeUnits::Milliseconds>() };
                TRY_AND_EXIT_ON_ITK_EXCEPTION(forwardProjection->Update());
            }

            image = forwardProjection->GetOutput();
        }
        else
        {
            {
                auto Timer{ HighPrecisionTimer<TimeUnits::Milliseconds>() };
                TRY_AND_EXIT_ON_ITK_EXCEPTION(resampleImageFilter->Update());
            }
            image = resampleImageFilter->GetOutput();
        }

        WriteImageFile(image, outputfilename, reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0]);

    }

    return EXIT_SUCCESS;
}

void WriteImageFile(OutputImageType::Pointer image, std::string & outputname, unsigned int size)
{
    outputname.append("_");
    outputname.append(std::to_string(size));
    outputname.append(".nrrd");
    using WriterType = itk::ImageFileWriter<  OutputImageType >;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputname);
    writer->SetInput(image);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update())
}
