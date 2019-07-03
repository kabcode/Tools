#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTranslationTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"

using PixelType = double;
static constexpr unsigned int Dim = 2;
using ImageType = itk::Image<PixelType, Dim>;
using ImageReaderType = itk::ImageFileReader<ImageType>;
using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType>;
using TranslationTransformType = itk::TranslationTransform<double, Dim>;

using MSMetricType = itk::MeanSquaresImageToImageMetric<ImageType, ImageType>;

template<int dimension>
void AssertDepencency()
{
	static_assert(dimension < 3, "function 2 not available");
	std::cout << "Here" << std::endl;
}

void AssertTest(int dimension)
{
	std::cout << dimension << std::endl;
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

	auto MSMetric = MSMetricType::New();
	MSMetric->SetInterpolator(Interpolator);
	MSMetric->SetTransform(TranslationTransform);
	MSMetric->SetFixedImage(FixedFileReader->GetOutput());
	MSMetric->SetFixedImageRegion(FixedFileReader->GetOutput()->GetLargestPossibleRegion());
	MSMetric->SetMovingImage(MovingFileReader->GetOutput());
	

	TranslationTransformType::ParametersType params(TranslationTransform->GetNumberOfParameters());
	params.Fill(0.0);

	MSMetric->Initialize();

	std::string OutputFileName("");
	if (argc < 4)
	{
		OutputFileName = "MetricMap.txt";
	}
	else
	{
		OutputFileName = argv[3];
	}

	std::ofstream outputFile(OutputFileName, std::ios::trunc);

	auto stepsize(1);

	if (!outputFile.is_open()) {
		std::cout << "File could not be opened." << std::endl;
	}
	else {
		for (auto x = -50.0; x <= 50.0; x += stepsize)
		{
			params(0) = x;
			for (auto y = -50.0; y <= 50.0; y += stepsize)
			{
				params(1) = y;
				try
				{
					auto value = MSMetric->GetValue(params);
					outputFile << x << "   " << y << "   " << value << std::endl;
					//std::cout << x << "\t" << y << " \t" << value << std::endl;
				}
				catch(...)
				{
					outputFile << x << "   " << y << "   " << "NaN" << std::endl;
					//std::cout << x << "\t" << y << " \t" << "Too many outside samples" << std::endl;
				}
				
				
			}
		}
	}

	outputFile.close();

	return EXIT_SUCCESS;


}