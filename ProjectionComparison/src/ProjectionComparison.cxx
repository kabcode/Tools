#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkJosephForwardProjectionImageFilter.h"
#include "rtkCudaForwardProjectionImageFilter.h"

#include "rtkConstantImageSource.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "HighPrecisionTimer.h"

using OutputPixelType = float;
constexpr unsigned int Dimension = 3;
using OutputImageType = itk::CudaImage< OutputPixelType, Dimension >;

void WriteImageFile(const OutputImageType::Pointer& image, std::string& outputname, unsigned int size, std::string additionalInformation = "");

int main(int argc, char * argv[])
{

    auto geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
    geometryReader->SetFilename(argv[2]);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation())

        // Create a stack of empty projection images
        using ConstantImageSourceType = rtk::ConstantImageSource< OutputImageType >;
    ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
    OutputImageType::SizeType OutputSize{ 2048,2048,1 };
    constantImageSource->SetSize(OutputSize);
    OutputImageType::SpacingType spacing;
    spacing.Fill(1);
    spacing[0] = 512. / OutputSize[0];
    spacing[1] = 512. / OutputSize[1];
    constantImageSource->SetSpacing(spacing);
    OutputImageType::PointType origin;
    origin[0] = -(OutputSize[0] * spacing[0]/2.0);
    origin[1] = -(OutputSize[1] * spacing[1] /2.0);
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

    rtk::ForwardProjectionImageFilter<OutputImageType>::Pointer forwardProjection;
    std::string outputfilename{ "" };

    /*
    forwardProjection = rtk::JosephForwardProjectionImageFilter<OutputImageType, OutputImageType>::New();
    outputfilename.append("Joseph");
    */
    
    forwardProjection = rtk::CudaForwardProjectionImageFilter<OutputImageType, OutputImageType>::New();
    dynamic_cast<rtk::CudaForwardProjectionImageFilter<OutputImageType, OutputImageType>*>(forwardProjection.GetPointer())->SetStepSize(1);
    outputfilename.append("Cuda");
    
    forwardProjection->SetInput(constantImageSource->GetOutput());
    forwardProjection->SetInput(1, reader->GetOutput());
    forwardProjection->SetGeometry(geometryReader->GetOutputObject());

    std::vector<long long> durations(105);
    for (unsigned int run = 0; run < durations.size(); ++run)
    {
        {
            auto Timer{ HighPrecisionTimer<TimeUnits::Microseconds, false>(&durations[run]) };
            TRY_AND_EXIT_ON_ITK_EXCEPTION(forwardProjection->Update());
        }
        forwardProjection->Modified();
    }
    const auto image = forwardProjection->GetOutput();
    WriteImageFile(image, outputfilename, reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0], "2048");


    durations.erase(durations.begin(), durations.begin() + 5);

    std::ofstream outstream{ "cuda_2048.txt", std::ios::out | std::ios::trunc };

    if (outstream.is_open())
    {
        for each (auto& datapoint in durations)
        {
            outstream << datapoint << std::endl;
        }
        outstream.close();
    }

    return EXIT_SUCCESS;
}

void WriteImageFile(const OutputImageType::Pointer& image, std::string & outputname, unsigned int size, std::string additionalInformation)
{
    outputname.append("_");
    outputname.append(std::to_string(size));
    if (!additionalInformation.empty())
    {
        outputname.append("_");
        outputname.append(additionalInformation);
    }

    outputname.append(".nrrd");
    using WriterType = itk::ImageFileWriter<  OutputImageType >;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputname);
    writer->SetInput(image);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update())
}
