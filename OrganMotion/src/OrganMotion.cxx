// STL includes
#include <filesystem>

// ITK includes
#include "itkImageFileReader.h"
#include "itkImageMomentsCalculator.h"

using PixelType = float;
const unsigned int Dim3D = 3;
using InputImageType = itk::Image<PixelType, Dim3D>;

using ImageReaderType = itk::ImageFileReader<InputImageType>;
using MomentCalculatorType = itk::ImageMomentsCalculator<InputImageType>;


int main(int argc, char* argv[])
{

    for (auto i = 1; i < argc; ++i)
    {
        auto ImageReader = ImageReaderType::New();
        ImageReader->SetFileName(argv[i]);

        auto Image = ImageReader->GetOutput();
        Image->Update();

        auto MomentCalculator = MomentCalculatorType::New();
        MomentCalculator->SetImage(Image);
        MomentCalculator->Compute();

        std::cout << i << ": " << MomentCalculator->GetCenterOfGravity() << std::endl;
    }
	


    

	return EXIT_SUCCESS;
	
}
