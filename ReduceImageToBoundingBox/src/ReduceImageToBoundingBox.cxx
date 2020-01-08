#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkLabelGeometryImageFilter.h"

using ImageType = itk::Image<float, 3>;


int main(int argc, char * argv[])
{
    using ImageReaderType = itk::ImageFileReader<ImageType>;
    ImageReaderType::Pointer labelReader = ImageReaderType::New();
    labelReader->SetFileName(argv[1]);
    labelReader->Update();

    auto labelImage = labelReader->GetOutput();

    auto intensityImage = ImageType::New();
    if(argc == 3)
    {
        ImageReaderType::Pointer intensityReader = ImageReaderType::New();
        intensityReader->SetFileName(argv[2]);
        intensityReader->Update();
        intensityImage = intensityReader->GetOutput();
    }
    
    using LabelGeometryImageFilterType = itk::LabelGeometryImageFilter<ImageType>;
    LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
    labelGeometryImageFilter->SetInput(labelImage);

    // These generate optional outputs.
    labelGeometryImageFilter->CalculatePixelIndicesOn();
    labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();
    labelGeometryImageFilter->CalculateOrientedLabelRegionsOn();
    
    labelGeometryImageFilter->Update();
    LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();

    LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
    std::cout << "Number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;
    std::cout << std::endl;

    for (allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++)
    {
        const LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
        std::cout << "\tLabel: " << static_cast<int>(labelValue) << std::endl;
        std::cout << "\tVolume: " << labelGeometryImageFilter->GetVolume(labelValue) << std::endl;
        std::cout << "\tCentroid: " << labelGeometryImageFilter->GetCentroid(labelValue) << std::endl;
        std::cout << "\tWeighted Centroid: " << labelGeometryImageFilter->GetWeightedCentroid(labelValue) << std::endl;
        std::cout << "\tAxes Length: " << labelGeometryImageFilter->GetAxesLength(labelValue) << std::endl;
        std::cout << "\tMajorAxisLength: " << labelGeometryImageFilter->GetMajorAxisLength(labelValue) << std::endl;
        std::cout << "\tMinorAxisLength: " << labelGeometryImageFilter->GetMinorAxisLength(labelValue) << std::endl;
        std::cout << "\tEccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue) << std::endl;
        std::cout << "\tElongation: " << labelGeometryImageFilter->GetElongation(labelValue) << std::endl;
        std::cout << "\tOrientation: " << labelGeometryImageFilter->GetOrientation(labelValue) << std::endl;
        std::cout << "\tBounding box: " << labelGeometryImageFilter->GetBoundingBox(labelValue) << std::endl;

        std::cout << std::endl << std::endl;
    }

    return EXIT_SUCCESS;
}
