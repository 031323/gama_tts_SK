/***************************************************************************
 *  Copyright 1991, 1992, 1993, 1994, 1995, 1996, 2001, 2002               *
 *    David R. Hill, Leonard Manzara, Craig Schock                         *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 ***************************************************************************/
// 2014-09
// This file was copied from Gnuspeech and modified by Marcelo Y. Matuda.

#ifndef VTM_VOCALTRACTMODEL0_H_
#define VTM_VOCALTRACTMODEL0_H_

#include <algorithm> /* max, min */
#include <istream>
#include <memory>
#include <vector>

#include "BandpassFilter.h"
#include "MovingAverageFilter.h"
#include "NoiseFilter.h"
#include "NoiseSource.h"
#include "RadiationFilter.h"
#include "ReflectionFilter.h"
#include "SampleRateConverter.h"
#include "Throat.h"
#include "VocalTractModelParameterValue.h"
#include "WavetableGlottalSource.h"

#define GS_VTM_TUBE_MIN_RADIUS (0.001)



namespace GS {
namespace VTM {

class VocalTractModel0 {
public:
	enum { /*  OROPHARYNX REGIONS  */
		R1 = 0, /*  S1  */
		R2 = 1, /*  S2  */
		R3 = 2, /*  S3  */
		R4 = 3, /*  S4 & S5  */
		R5 = 4, /*  S6 & S7  */
		R6 = 5, /*  S8  */
		R7 = 6, /*  S9  */
		R8 = 7, /*  S10  */
		TOTAL_REGIONS = 8
	};
	enum { /*  NASAL TRACT SECTIONS  */
		N1 = 0,
		N2 = 1,
		N3 = 2,
		N4 = 3,
		N5 = 4,
		N6 = 5,
		TOTAL_NASAL_SECTIONS = 6
	};
	enum ParameterIndex {
		PARAM_GLOT_PITCH = 0,
		PARAM_GLOT_VOL   = 1,
		PARAM_ASP_VOL    = 2,
		PARAM_FRIC_VOL   = 3,
		PARAM_FRIC_POS   = 4,
		PARAM_FRIC_CF    = 5,
		PARAM_FRIC_BW    = 6,
		PARAM_R1         = 7,
		PARAM_R2         = 8,
		PARAM_R3         = 9,
		PARAM_R4         = 10,
		PARAM_R5         = 11,
		PARAM_R6         = 12,
		PARAM_R7         = 13,
		PARAM_R8         = 14,
		PARAM_VELUM      = 15
	};
	enum Waveform {
		GLOTTAL_SOURCE_PULSE = 0,
		GLOTTAL_SOURCE_SINE = 1
	};

	struct Configuration {
		float  outputRate;                  /*  output sample rate (22.05, 44.1)  */
		float  controlRate;                 /*  1.0-1000.0 input tables/second (Hz)  */

		double volume;                      /*  master volume (0 - 60 dB)  */
		int    channels;                    /*  # of sound output channels (1, 2)  */
		double balance;                     /*  stereo balance (-1 to +1)  */

		int    waveform;                    /*  GS waveform type (0=PULSE, 1=SINE  */
		double tp;                          /*  % glottal pulse rise time  */
		double tnMin;                       /*  % glottal pulse fall time minimum  */
		double tnMax;                       /*  % glottal pulse fall time maximum  */
		double breathiness;                 /*  % glottal source breathiness  */

		double length;                      /*  nominal tube length (10 - 20 cm)  */
		double temperature;                 /*  tube temperature (25 - 40 C)  */
		double lossFactor;                  /*  junction loss factor in (0 - 5 %)  */

		double apertureRadius;              /*  aperture scl. radius (3.05 - 12 cm)  */
		double mouthCoef;                   /*  mouth aperture coefficient  */
		double noseCoef;                    /*  nose aperture coefficient  */

		// Set noseRadius[N1] to 0.0, because it is not used.
		double noseRadius[TOTAL_NASAL_SECTIONS]; /*  fixed nose radii (0 - 3 cm)  */

		double throatCutoff;                /*  throat lp cutoff (50 - nyquist Hz)  */
		double throatVol;                   /*  throat volume (0 - 48 dB) */

		int    modulation;                  /*  pulse mod. of noise (0=OFF, 1=ON)  */
		double mixOffset;                   /*  noise crossmix offset (30 - 60 dB)  */
	};

	VocalTractModel0();
	~VocalTractModel0();

	void synthesizeToFile(std::istream& inputStream, const char* outputFile);
	void synthesizeToBuffer(std::istream& inputStream, std::vector<float>& outputBuffer);

	void loadConfiguration(const Configuration& config) {
		config_ = config;
	}

	void initializeSynthesizer();
	void initializeInputFilters(double period);
	void loadSingleInput(const VocalTractModelParameterValue pv);
	void synthesizeForInputSequence();
	void synthesizeForSingleInput(int numIterations);

	std::vector<float>& outputData() { return outputData_; }
	std::size_t outputDataPos() const { return outputDataPos_; }
	void setOutputDataPos(std::size_t pos) { outputDataPos_ = pos; }
	void resetOutputData() {
		outputData_.clear();
		outputDataPos_ = 0;
	}
	double maximumOutputSampleValue() const { return srConv_->maximumSampleValue(); }
	float outputRate() const { return config_.outputRate; }
	unsigned int numChannels() const { return config_.channels; }
private:
	enum {
		VELUM = N1
	};
	enum { /*  OROPHARYNX SCATTERING JUNCTION COEFFICIENTS (BETWEEN EACH REGION)  */
		C1 = R1, /*  R1-R2 (S1-S2)  */
		C2 = R2, /*  R2-R3 (S2-S3)  */
		C3 = R3, /*  R3-R4 (S3-S4)  */
		C4 = R4, /*  R4-R5 (S5-S6)  */
		C5 = R5, /*  R5-R6 (S7-S8)  */
		C6 = R6, /*  R6-R7 (S8-S9)  */
		C7 = R7, /*  R7-R8 (S9-S10)  */
		C8 = R8, /*  R8-AIR (S10-AIR)  */
		TOTAL_COEFFICIENTS = TOTAL_REGIONS
	};
	enum { /*  OROPHARYNX SECTIONS  */
		S1  = 0, /*  R1  */
		S2  = 1, /*  R2  */
		S3  = 2, /*  R3  */
		S4  = 3, /*  R4  */
		S5  = 4, /*  R4  */
		S6  = 5, /*  R5  */
		S7  = 6, /*  R5  */
		S8  = 7, /*  R6  */
		S9  = 8, /*  R7  */
		S10 = 9, /*  R8  */
		TOTAL_SECTIONS = 10
	};
	enum { /*  NASAL TRACT COEFFICIENTS  */
		NC1 = N1, /*  N1-N2  */
		NC2 = N2, /*  N2-N3  */
		NC3 = N3, /*  N3-N4  */
		NC4 = N4, /*  N4-N5  */
		NC5 = N5, /*  N5-N6  */
		NC6 = N6, /*  N6-AIR  */
		TOTAL_NASAL_COEFFICIENTS = TOTAL_NASAL_SECTIONS
	};
	enum { /*  THREE-WAY JUNCTION ALPHA COEFFICIENTS  */
		LEFT  = 0,
		RIGHT = 1,
		UPPER = 2,
		TOTAL_ALPHA_COEFFICIENTS = 3
	};
	enum { /*  FRICATION INJECTION COEFFICIENTS  */
		FC1 = 0, /*  S3  */
		FC2 = 1, /*  S4  */
		FC3 = 2, /*  S5  */
		FC4 = 3, /*  S6  */
		FC5 = 4, /*  S7  */
		FC6 = 5, /*  S8  */
		FC7 = 6, /*  S9  */
		FC8 = 7, /*  S10  */
		TOTAL_FRIC_COEFFICIENTS = 8
	};

	struct InputData {
		double glotPitch;
		double glotVol;
		double aspVol;
		double fricVol;
		double fricPos;
		double fricCF;
		double fricBW;
		double radius[TOTAL_REGIONS];
		double velum;
	};

	/*  VARIABLES FOR INTERPOLATION  */
	struct CurrentData {
		double glotPitch;
		double glotPitchDelta;
		double glotVol;
		double glotVolDelta;
		double aspVol;
		double aspVolDelta;
		double fricVol;
		double fricVolDelta;
		double fricPos;
		double fricPosDelta;
		double fricCF;
		double fricCFDelta;
		double fricBW;
		double fricBWDelta;
		double radius[TOTAL_REGIONS];
		double radiusDelta[TOTAL_REGIONS];
		double velum;
		double velumDelta;
	};

	struct InputFilters {
		MovingAverageFilter<double> glotPitchFilter;
		MovingAverageFilter<double> glotVolFilter;
		MovingAverageFilter<double> aspVolFilter;
		MovingAverageFilter<double> fricVolFilter;
		MovingAverageFilter<double> fricPosFilter;
		MovingAverageFilter<double> fricCFFilter;
		MovingAverageFilter<double> fricBWFilter;
		MovingAverageFilter<double> radius0Filter;
		MovingAverageFilter<double> radius1Filter;
		MovingAverageFilter<double> radius2Filter;
		MovingAverageFilter<double> radius3Filter;
		MovingAverageFilter<double> radius4Filter;
		MovingAverageFilter<double> radius5Filter;
		MovingAverageFilter<double> radius6Filter;
		MovingAverageFilter<double> radius7Filter;
		MovingAverageFilter<double> velumFilter;
		InputFilters(double sampleRate, double period)
			: glotPitchFilter(sampleRate, period)
			, glotVolFilter(sampleRate, period)
			, aspVolFilter(sampleRate, period)
			, fricVolFilter(sampleRate, period)
			, fricPosFilter(sampleRate, period)
			, fricCFFilter(sampleRate, period)
			, fricBWFilter(sampleRate, period)
			, radius0Filter(sampleRate, period)
			, radius1Filter(sampleRate, period)
			, radius2Filter(sampleRate, period)
			, radius3Filter(sampleRate, period)
			, radius4Filter(sampleRate, period)
			, radius5Filter(sampleRate, period)
			, radius6Filter(sampleRate, period)
			, radius7Filter(sampleRate, period)
			, velumFilter(sampleRate, period) {}
		void reset() {
			glotPitchFilter.reset();
			glotVolFilter.reset();
			aspVolFilter.reset();
			fricVolFilter.reset();
			fricPosFilter.reset();
			fricCFFilter.reset();
			fricBWFilter.reset();
			radius0Filter.reset();
			radius1Filter.reset();
			radius2Filter.reset();
			radius3Filter.reset();
			radius4Filter.reset();
			radius5Filter.reset();
			radius6Filter.reset();
			radius7Filter.reset();
			velumFilter.reset();
		}
	};

	VocalTractModel0(const VocalTractModel0&) = delete;
	VocalTractModel0& operator=(const VocalTractModel0&) = delete;

	void reset();
	void calculateTubeCoefficients();
	void initializeNasalCavity();
	void parseInputStream(std::istream& in);
	void sampleRateInterpolation();
	void setControlRateParameters(int pos);
	void setFricationTaps();
	double vocalTract(double input, double frication);
	void writeOutputToFile(const char* outputFile);
	void writeOutputToBuffer(std::vector<float>& outputBuffer);
	void synthesize();
	float calculateMonoScale();
	void calculateStereoScale(float& leftScale, float& rightScale);

	template<typename T> T readParameterFromInputStream(std::istream& in, std::string& line, const char* paramName);

	Configuration config_;

	/*  DERIVED VALUES  */
	int    controlPeriod_;
	int    sampleRate_;
	double actualTubeLength_;            /*  actual length in cm  */

	/*  MEMORY FOR TUBE AND TUBE COEFFICIENTS  */
	double oropharynx_[TOTAL_SECTIONS][2][2];
	double oropharynxCoeff_[TOTAL_COEFFICIENTS];

	double nasal_[TOTAL_NASAL_SECTIONS][2][2];
	double nasalCoeff_[TOTAL_NASAL_COEFFICIENTS];

	double alpha_[TOTAL_ALPHA_COEFFICIENTS];
	int currentPtr_;
	int prevPtr_;

	/*  MEMORY FOR FRICATION TAPS  */
	double fricationTap_[TOTAL_FRIC_COEFFICIENTS];

	double dampingFactor_;               /*  calculated damping factor  */
	double crossmixFactor_;              /*  calculated crossmix factor  */
	double breathinessFactor_;

	double prevGlotAmplitude_;

	std::vector<std::unique_ptr<InputData>> inputData_;
	CurrentData currentData_;
	InputData singleInput_;
	std::size_t outputDataPos_;
	std::vector<float> outputData_;
	std::unique_ptr<SampleRateConverter> srConv_;
	std::unique_ptr<RadiationFilter> mouthRadiationFilter_;
	std::unique_ptr<ReflectionFilter> mouthReflectionFilter_;
	std::unique_ptr<RadiationFilter> nasalRadiationFilter_;
	std::unique_ptr<ReflectionFilter> nasalReflectionFilter_;
	std::unique_ptr<Throat> throat_;
	std::unique_ptr<WavetableGlottalSource> glottalSource_;
	std::unique_ptr<BandpassFilter> bandpassFilter_;
	std::unique_ptr<NoiseFilter> noiseFilter_;
	std::unique_ptr<NoiseSource> noiseSource_;
	std::unique_ptr<InputFilters> inputFilters_;
};

} /* namespace VTM */
} /* namespace GS */

#endif /* VTM_VOCALTRACTMODEL0_H_ */
