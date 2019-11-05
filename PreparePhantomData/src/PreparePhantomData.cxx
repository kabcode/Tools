#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"


using ImageType = itk::Image<float, 3>;

int main(int argc, char* argv[])
{
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

    auto AddImageFilter = itk::AddImageFilter<ImageType>::New();
    AddImageFilter->SetInput1(ImageReader->GetOutput());
    AddImageFilter->SetConstant2(1000);

    auto MultiplyImageFilter = itk::MultiplyImageFilter<ImageType>::New();
    MultiplyImageFilter->SetInput1(AddImageFilter->GetOutput());
    MultiplyImageFilter->SetConstant2(100);

    auto MinimumMaximumImageFilter = itk::MinimumMaximumImageFilter<ImageType>::New();
    MinimumMaximumImageFilter->SetInput(MultiplyImageFilter->GetOutput());
    MinimumMaximumImageFilter->Update();

    auto min = MinimumMaximumImageFilter->GetMinimum();
    auto max = MinimumMaximumImageFilter->GetMaximum();



}
