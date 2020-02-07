#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTranslationTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkAdvGradientDifferenceImageToImageMetric.h"

using PixelType = double;
static constexpr unsigned int Dim = 2;
using ImageType = itk::Image<PixelType, Dim>;
using ImageReaderType = itk::ImageFileReader<ImageType>;
using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType>;
using TranslationTransformType = itk::TranslationTransform<double, Dim>;

using AGDMetricType = itk::AdvGradientDifferenceImageToImageMetric<ImageType, ImageType>;

template<int dimension>
void AssertDepencency()
{
	static_assert(dimension < 3, "function 2 not available");
}

void AssertTest(int dimension)
{
	AssertDepencency<2>();
}

void
PrintUsage(std::string programme)
{
	std::cout << "USAGE:\n";
	std::cout << programme << " FixedImageFile MovingImageFile \n";
	std::cout << "Output:" << " MetricMap.txt with 3 columns (x,y,metric value) \n";
}

inline int computeSteps(std::pair<int,int> range, int stepsize)
{
    return std::abs(range.second - range.first)/(double)stepsize + 1;
}

int main(int argc, char* argv[])
{
	if (argc < 3 || argc > 4)
	{
		PrintUsage(argv[0]);
		return EXIT_FAILURE;
	}

	AssertTest(2);

	auto FixedFileReader = ImageReaderType::New();
	FixedFileReader->SetFileName(argv[1]);

	auto MovingFileReader = ImageReaderType::New();
	MovingFileReader->SetFileName(argv[2]);

	try
	{
		FixedFileReader->Update();
		MovingFileReader->Update();
	}
	catch(itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
		return EXIT_FAILURE;
	}


	auto TranslationTransform = TranslationTransformType::New();
	auto Interpolator = InterpolatorType::New();
	Interpolator->SetInputImage(FixedFileReader->GetOutput());

	auto AGDMetric = AGDMetricType::New();
    AGDMetric->SetInterpolator(Interpolator);
    AGDMetric->SetTransform(TranslationTransform);
    AGDMetric->SetFixedImage(FixedFileReader->GetOutput());
    AGDMetric->SetFixedImageRegion(FixedFileReader->GetOutput()->GetLargestPossibleRegion());
    AGDMetric->SetMovingImage(MovingFileReader->GetOutput());
	

	TranslationTransformType::ParametersType params(TranslationTransform->GetNumberOfParameters());
	params.Fill(0.0);

    AGDMetric->Initialize();

	std::string OutputFileName;
	if (argc < 4)
	{
		OutputFileName = "MetricMap.txt";
	}
	else
	{
		OutputFileName = argv[3];
	}

	std::ofstream outputFile(OutputFileName, std::ios::trunc);

	const auto stepsize(0.5);
    const std::pair<int, int> x_range{ -60,60 };
    const std::pair<int, int> y_range{ -400,400 };

    std::pair<int, int> steps{computeSteps(x_range, stepsize), computeSteps(y_range, stepsize)};

    outputFile << "0 " << steps.first << " 1 " << steps.second << " " << stepsize << std::endl;

	if (!outputFile.is_open()) {
		std::cout << "File could not be opened." << std::endl;
	}
	else {
		for (auto x = x_range.first; x <= x_range.second; x += stepsize)
		{
			params(0) = x;
			for (auto y = y_range.first; y <= y_range.second; y += stepsize)
			{
				params(1) = y;
				try
				{
					const auto value = AGDMetric->GetValue(params);
                    outputFile <<"[" << x << ", " << y << "], " << value << std::endl;
					
				}
				catch(...)
				{
                    outputFile << "[" << x << ", " << y << "], " << "Nan" << std::endl;
				}
				
				
			}
		}
	}

	outputFile.close();

	return EXIT_SUCCESS;


}