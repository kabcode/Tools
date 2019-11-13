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

    std::cout << AGDMetric->GetValue(params) << std::endl;

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

	const auto stepsize(5);

	if (!outputFile.is_open()) {
		std::cout << "File could not be opened." << std::endl;
	}
	else {
		for (auto x = -20; x <= 20; x += stepsize)
		{
			params(0) = x;
			for (auto y = -20; y <= 20; y += stepsize)
			{
				params(1) = y;
				try
				{
					const auto value = AGDMetric->GetValue(params);
                    outputFile << value << " ";
					//std::cout << x << "\t" << y << " \t" << value << std::endl;
				}
				catch(...)
				{
					outputFile << "NaN" << " ";
					//std::cout << x << "\t" << y << " \t" << "Too many outside samples" << std::endl;
				}
				
				
			}
            outputFile << std::endl;
		}
	}

	outputFile.close();

	return EXIT_SUCCESS;


}