#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkRegionOfInterestImageFilter.h"


using ImageType = itk::Image<float, 3>;
using ImageToHistogrammFilterType = itk::Statistics::ImageToHistogramFilter<ImageType>;

int main(int argc, char* argv[])
{
    std::cout << "Load image:" << argv[1] << std::endl;
    auto ImageReader = itk::ImageFileReader<ImageType>::New();
    ImageReader->SetFileName(argv[1]);

    try
    {
        ImageReader->Update();
    }
    catch (itk::ExceptionObject &EO)
    {
        std::cout << "Read failed for: " << argv << std::endl;
        std::cout << EO << std::endl;
    }

    std::cout << "AddImageFilter" << std::endl;
    auto AddImageFilter = itk::AddImageFilter<ImageType>::New();
    AddImageFilter->SetInput1(ImageReader->GetOutput());
    AddImageFilter->SetConstant2(1000);

    std::cout << "MultiplyImageFilter" << std::endl;
    auto MultiplyImageFilter = itk::MultiplyImageFilter<ImageType>::New();
    MultiplyImageFilter->SetInput1(AddImageFilter->GetOutput());
    MultiplyImageFilter->SetConstant2(100);

    std::cout << "MinMaxImageFilter" << std::endl;
    auto MinimumMaximumImageFilter = itk::MinimumMaximumImageFilter<ImageType>::New();
    MinimumMaximumImageFilter->SetInput(MultiplyImageFilter->GetOutput());
    MinimumMaximumImageFilter->Update();

    const auto min = MinimumMaximumImageFilter->GetMinimum();
    const auto max = MinimumMaximumImageFilter->GetMaximum();
    std::cout << "\tMin:" << min << std::endl;
    std::cout << "\tMax:" << max << std::endl;

    std::cout << "ImagetoHistogrammFilter" << std::endl;
    auto ImagetoHistogrammFilter = ImageToHistogrammFilterType::New();
    ImageToHistogrammFilterType::HistogramSizeType HistogrammArray(1);
    HistogrammArray[0] = max - min;
    ImagetoHistogrammFilter->SetHistogramSize(HistogrammArray);
    ImagetoHistogrammFilter->SetMarginalScale(1);
    ImagetoHistogrammFilter->SetAutoMinimumMaximum(true);

    ImagetoHistogrammFilter->SetInput(MultiplyImageFilter->GetOutput());

    try
    {
        ImagetoHistogrammFilter->Update();
    }
    catch (itk::ExceptionObject &EO)
    {
        EO.Print(std::cout);
        return EXIT_FAILURE;
    }

    std::cout << "Histogram of image" << std::endl;
    const ImageToHistogrammFilterType::HistogramType * histogram = ImagetoHistogrammFilter->GetOutput();
    for (unsigned long long bin = 0; bin < HistogrammArray[0] ; ++bin)
    {
        if(histogram->GetFrequency(bin, 0) > 0)
        {
            auto ThresholdImageFilter = itk::ThresholdImageFilter<ImageType>::New();
            ThresholdImageFilter->SetInput(MultiplyImageFilter->GetOutput());
            ThresholdImageFilter->SetOutsideValue(0);
            ThresholdImageFilter->SetLower(histogram->GetFrequency(bin, 0) - 1);
            ThresholdImageFilter->SetLower(histogram->GetFrequency(bin, 0) + 1);
        } 
    }

    auto LabelImageToShapeLabelMapFilter = itk::LabelImageToShapeLabelMapFilter<ImageType>::New();
    LabelImageToShapeLabelMapFilter->SetInput(MultiplyImageFilter->GetOutput());
    try
    {
        LabelImageToShapeLabelMapFilter->Update();
    }
    catch (itk::ExceptionObject &EO)
    {
        EO.Print(std::cout);
        return EXIT_FAILURE;
    }

    const auto LabelMapPointer = LabelImageToShapeLabelMapFilter->GetOutput();
    std::cout << "Size: " << LabelMapPointer->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "Number: " << LabelMapPointer->GetNumberOfLabelObjects() << std::endl;

    for (unsigned long long i = 0; i < LabelMapPointer->GetNumberOfLabelObjects(); ++i)
    {
        const auto LabelObject = LabelMapPointer->GetNthLabelObject(i);
        auto MaskedImage = ImageType::New();
        auto region = ImageReader->GetOutput()->GetLargestPossibleRegion();
        MaskedImage->SetRegions(region);
        MaskedImage->Allocate();
        MaskedImage->FillBuffer(0);

        MaskedImage->SetSpacing(ImageReader->GetOutput()->GetSpacing());
        MaskedImage->SetDirection(ImageReader->GetOutput()->GetDirection());
        MaskedImage->SetOrigin(ImageReader->GetOutput()->GetOrigin());

        itk::ImageRegionIterator<ImageType> MaskIterator(MaskedImage, MaskedImage->GetLargestPossibleRegion());
        itk::ImageRegionConstIterator<ImageType> LabelIterator(ImageReader->GetOutput(), ImageReader->GetOutput()->GetLargestPossibleRegion());

        /* Copy the label into a new image */
        while (!MaskIterator.IsAtEnd())
        {
            if (LabelIterator.Get() == i)
            {
                MaskIterator.Set(LabelIterator.Get());
            }
            ++MaskIterator;
            ++LabelIterator;
        }


        /* Crop masked intensity images to the bounding box */
        auto RegionOfInterestFilter = itk::RegionOfInterestImageFilter<ImageType, ImageType>::New();
        RegionOfInterestFilter->SetInput(MaskedImage);
        RegionOfInterestFilter->SetRegionOfInterest(LabelObject->GetBoundingBox());

        std::string FileName("Label_");
        FileName.append(std::to_string(i));
        FileName.append(".nrrd");
        auto ImageWriter = itk::ImageFileWriter<ImageType>::New();
        ImageWriter->SetInput(RegionOfInterestFilter->GetOutput());
        ImageWriter->SetFileName(FileName);

        try
        {
            ImageWriter->Update();
        }
        catch (itk::ExceptionObject &EO)
        {
            EO.Print(std::cout);
            return EXIT_FAILURE;
        }
    }

    std::cout << "WriteImage" << std::endl;
    auto ImageFileWriter = itk::ImageFileWriter<ImageType>::New();
    ImageFileWriter->SetInput(MultiplyImageFilter->GetOutput());
    ImageFileWriter->SetFileName("Output.nrrd");

    try
    {
        ImageFileWriter->Update();
    }
    catch (itk::ExceptionObject &EO)
    {
        EO.Print(std::cout);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
