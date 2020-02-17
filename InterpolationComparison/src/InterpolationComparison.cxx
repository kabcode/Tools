#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "HighPrecisionTimer.h"
#include "itkSubtractImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkMacro.h"
#include "rtkMacro.h"

#include "itkResampleImageFilter.h"
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkGaussianInterpolateImageFunction.h>

using OutputPixelType = double;
constexpr unsigned int Dimension = 3;
using ImageType = itk::Image< OutputPixelType, Dimension >;
using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType, double>;
using NearestNeighborInterpolationType = itk::NearestNeighborInterpolateImageFunction<ImageType, double>;
using LinearInterpolationType = itk::LinearInterpolateImageFunction<ImageType, double>;
using BSlineInterpolationType = itk::BSplineInterpolateImageFunction<ImageType, double>;
using GaussianInterpolationType = itk::GaussianInterpolateImageFunction<ImageType, double>;
using ExtractImageFilterType = itk::ExtractImageFilter<ImageType, itk::Image<double, 2>>;

template<class TImage>
void WriteImageFile(const typename TImage::Pointer image, std::string outputname, std::string additionalInformation = "");

auto Extract2DImage(ImageType::Pointer image) ->itk::Image<double, 2>::Pointer;

int main(int argc, char * argv[])
{
    auto reader = itk::ImageFileReader<ImageType>::New();
    reader->SetFileName(argv[1]);

    auto InputImage = reader->GetOutput();
    TRY_AND_EXIT_ON_ITK_EXCEPTION(InputImage->Update());

    auto Transform = itk::MatrixOffsetTransformBase<double>::New();
    Transform->SetIdentity();

    auto ResampleImageFilter = ResampleImageFilterType::New();

    ImageType::SizeType ImageSize;
    ImageSize.Fill(1024);
    ImageSize[2] = 1;
    ImageType::SpacingType ImageSpacing;
    ImageSpacing.Fill(0.5);
    ImageSpacing[2] = 0.5;
    ImageType::PointType Origin;
    Origin.Fill(-255.75);
    ResampleImageFilter->SetInput(InputImage);
    ResampleImageFilter->SetSize(ImageSize);
    ResampleImageFilter->SetOutputDirection(InputImage->GetDirection());
    ResampleImageFilter->SetOutputOrigin(Origin);
    ResampleImageFilter->SetOutputSpacing(ImageSpacing);
    ResampleImageFilter->SetTransform(Transform);

    itk::InterpolateImageFunction<ImageType, double>::Pointer Interpolator = nullptr;

    for (unsigned int TYPE = 0; TYPE < 4; ++TYPE)
    {
        std::string OutputFilename{};
        switch (TYPE)
        {
        case 0:
            Interpolator = NearestNeighborInterpolationType::New();
            OutputFilename.append("NearestNeigbor");
            break;
        case 1:
            Interpolator = LinearInterpolationType::New();
            OutputFilename.append("Linear");
            break;
        case 2:
            Interpolator = BSlineInterpolationType::New();
            OutputFilename.append("BSPline");
            break;
        case 3:
            Interpolator = GaussianInterpolationType::New();
            OutputFilename.append("Gaussian");
            break;
        default:
            std::cout << "Value not handled - abort";
            return EXIT_FAILURE;
        }

        if (dynamic_cast<itk::GaussianInterpolateImageFunction<ImageType, double>*>(Interpolator.GetPointer()) != nullptr)
        {
            auto nInterpolator = dynamic_cast<itk::GaussianInterpolateImageFunction<ImageType, double>*>(Interpolator.GetPointer());
            nInterpolator->SetSigma(0.8);
            Interpolator = nInterpolator;
        }

        ResampleImageFilter->SetInterpolator(Interpolator);

        std::vector<long long> durations(105);
        for (long long& duration : durations)
        {
            {
                auto Timer{ HighPrecisionTimer<TimeUnits::Milliseconds, false>(&duration) };
                TRY_AND_EXIT_ON_ITK_EXCEPTION(ResampleImageFilter->Update());
            }
            ResampleImageFilter->Modified();
        }
        durations.erase(durations.begin(), durations.begin() + 5);

        std::ofstream outstream{ std::string{OutputFilename}.append(".txt"), std::ios::out | std::ios::trunc };

        if (outstream.is_open())
        {
            for each (auto& datapoint in durations)
            {
                outstream << datapoint << std::endl;
            }
            outstream.close();
        }

        WriteImageFile<ImageType>(ResampleImageFilter->GetOutput(), OutputFilename);

        auto GTReader = itk::ImageFileReader<ImageType>::New();
        GTReader->SetFileName(argv[2]);

        auto SubtractionImageFilter = itk::SubtractImageFilter<ImageType>::New();
        SubtractionImageFilter->SetInput1(ResampleImageFilter->GetOutput());
        SubtractionImageFilter->SetInput2(GTReader->GetOutput());
        auto AbsImageFilter = itk::AbsImageFilter<ImageType, ImageType>::New();
        AbsImageFilter->SetInput(SubtractionImageFilter->GetOutput());
        TRY_AND_EXIT_ON_ITK_EXCEPTION(AbsImageFilter->Update());

        WriteImageFile<itk::Image<double, 2>>(Extract2DImage(AbsImageFilter->GetOutput()), OutputFilename, "Sub");
        WriteImageFile<itk::Image<double, 2>>(Extract2DImage(GTReader->GetOutput()), "GroundTruth", "2D");
    }
    return EXIT_SUCCESS;
}

template<class TImage>
void WriteImageFile(typename TImage::Pointer image, std::string outputname, std::string additionalInformation)
{
    if (!additionalInformation.empty())
    {
        outputname.append("_");
        outputname.append(additionalInformation);
    }

    outputname.append(".nrrd");
    using WriterType = itk::ImageFileWriter<TImage>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputname);
    writer->SetInput(image);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update())
}

auto Extract2DImage(ImageType::Pointer image)-> itk::Image<double,2>::Pointer
{
    auto ExtractImageFilter = ExtractImageFilterType::New();
    auto region{ image->GetLargestPossibleRegion() };
    auto size = region.GetSize();
    size[2] = 0;
    region.SetSize(size);
    ExtractImageFilter->SetExtractionRegion(region);
    ExtractImageFilter->SetDirectionCollapseToIdentity();
    ExtractImageFilter->SetInput(image);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(ExtractImageFilter->Update());

    return ExtractImageFilter->GetOutput();
}