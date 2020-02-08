#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"

const unsigned int Dimension = 3;
using PixelType = float;
using ImageType = itk::Image<PixelType, Dimension>;


using ImageFileReaderType = itk::ImageFileReader<ImageType>;
using ImageFileWriterType = itk::ImageFileWriter<ImageType>;

inline auto CheckSpacing(ImageType::SpacingValueType in, ImageType::SpacingValueType sp) ->ImageType::SpacingValueType
{
    return (in == -1.f) ? sp : in;
}
inline auto ComputeOutputOrigin(ImageType::PointType or ,ImageType::SpacingType in, ImageType::SpacingType out) -> ImageType::PointType
{
    ImageType::PointType outputOrigin;
    for(unsigned int i = 0; i < or.Size(); ++i)
    {
        if(out[i] < in[i])
            outputOrigin[i] = or [i] + 0.5 * out[i] / in[i];
        
        if (out[i] > in[i])
            outputOrigin[i] = or [i] - 0.5 * out[i] / in[i];

        if (out[i] == in[i])
            outputOrigin[i] = or[i];
    }
    return outputOrigin;
}

inline auto ResizeDimension(ImageType::SizeType in_sz, ImageType::SpacingType in_sp, ImageType::SpacingType out_sp) -> ImageType::SizeType
{
    ImageType::SizeType out_sz;
    for(unsigned int i = 0; i < in_sz.size(); ++i)
    {
        out_sz[i] = in_sz[i] / out_sp[i] * in_sp[i];
    }
    return out_sz;
}

auto WriteImage(const ImageType::Pointer, const std::string&) -> void;

int main(int argc, char *argv[])
{
	if(argc < 6)
	{
		std::cout << argv[0] << " InputImageFileName -o OutputFilename -sp x y z" << std::endl;
		return EXIT_FAILURE;
	}

    auto ImageFileReader = ImageFileReaderType::New();
    ImageFileReader->SetFileName(argv[1]);
    auto ResamplingImage = ImageFileReader->GetOutput();
    try
    {
        ResamplingImage->Update();
    }
    catch (itk::ExceptionObject &EO)
    {
        EO.Print(std::cout);
        return EXIT_FAILURE;
    }

    std::string OutputFilename{ "" };
    auto outputSpacing = ImageFileReader->GetOutput()->GetSpacing();
	for (auto i = 1; i < argc; ++i)
	{
		if (std::string(argv[i]) == "-o")
		{
			++i;
			OutputFilename = argv[i];
		}
        if (std::string(argv[i]) == "-sp")
        {
            ++i;
            outputSpacing[0] = CheckSpacing(std::strtof(argv[i],nullptr), outputSpacing[0]);
            ++i;
            outputSpacing[1] = CheckSpacing(std::strtof(argv[i], nullptr), outputSpacing[1]);
            ++i;
            outputSpacing[2] = CheckSpacing(std::strtof(argv[i], nullptr), outputSpacing[2]);
        }
	}

    auto ResampleFilter = itk::ResampleImageFilter<ImageType, ImageType>::New();
    ResampleFilter->SetInput(ResamplingImage);

    const auto LinearInt = itk::LinearInterpolateImageFunction<ImageType>::New();
    ResampleFilter->SetInterpolator(LinearInt);
  
    ResampleFilter->SetOutputDirection(ImageFileReader->GetOutput()->GetDirection());
    const auto OutputOrigin = ComputeOutputOrigin(ImageFileReader->GetOutput()->GetOrigin(), ImageFileReader->GetOutput()->GetSpacing(), outputSpacing);
    ResampleFilter->SetOutputOrigin(OutputOrigin);
    const auto OutputSize = ResizeDimension(ResamplingImage->GetLargestPossibleRegion().GetSize(), ResamplingImage->GetSpacing(), outputSpacing);
    ResampleFilter->SetSize(OutputSize);
    ResampleFilter->SetOutputSpacing(outputSpacing);
    ResampleFilter->Update();

    WriteImage(ResampleFilter->GetOutput(), OutputFilename);
	
	return EXIT_SUCCESS;
}

void
WriteImage(const ImageType::Pointer image, const std::string& outfilename)
{
	auto ImageFileWriter = ImageFileWriterType::New();
	ImageFileWriter->SetInput(image);
	std::string OutputImageFilename(outfilename);
	ImageFileWriter->SetFileName(OutputImageFilename);

	try
	{
		ImageFileWriter->Update();
	}
	catch (itk::ExceptionObject &EO)
	{
		EO.Print(std::cout);
	}
}