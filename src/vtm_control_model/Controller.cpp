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

#include "Controller.h"

#include <cctype> /* isspace */
#include <cmath> /* rint */
#include <cstdio> /* printf */
#include <fstream>
#include <sstream>

#include "Exception.h"
#include "Log.h"
#include "vtm_plugin.h"
#include "VTMUtil.h"
#include "WAVEFileWriter.h"

#define VTM_CONFIG_FILE_NAME "/vtm.config"



namespace GS {
namespace VTMControlModel {

Controller::Controller(const char* configDirPath, Model& model,bool sks)
		: configDirPath_(configDirPath)
		, model_(model)
		, eventList_(configDirPath, model_)
		, vtmControlModelConfig_(configDirPath)
		, outputScale_(1.0)
		, sk(sks)

{
	// Load VTM configuration.
	std::ostringstream vtmConfigFilePath;
	vtmConfigFilePath << configDirPath << VTM_CONFIG_FILE_NAME;
	vtmConfigData_ = std::make_unique<ConfigurationData>(vtmConfigFilePath.str());
	vtmConfigData_->insert(*vtmControlModelConfig_.voiceData);

	// Get the vocal tract model instance.
	vtm_ = std::make_unique<VTM::VocalTractModelPlugin>(*vtmConfigData_);

	eventList_.setControlPeriod(vtmControlModelConfig_.controlPeriod);
}

void
Controller::initUtterance()
{
	int outputRate = vtmConfigData_->value<int>("output_rate");
	const float vtlOffset = vtmConfigData_->value<float>("vocal_tract_length_offset");
	const float vocalTractLength = vtmConfigData_->value<float>("vocal_tract_length");

	if (Log::debugEnabled) {
		printf("Tube Length = %f\n", vtlOffset + vocalTractLength);
		printf("Voice: %s\n", vtmControlModelConfig_.voiceName.c_str());
		printf("sampling Rate: %d\n", outputRate);
	}

	eventList_.setInitialPitch(vtmControlModelConfig_.initialPitch);
	eventList_.setMeanPitch(vtmControlModelConfig_.pitchOffset + vtmConfigData_->value<double>("reference_glottal_pitch"));
	eventList_.setGlobalTempo(vtmControlModelConfig_.tempo);
	eventList_.setUpDriftGenerator(vtmControlModelConfig_.driftDeviation, vtmControlModelConfig_.controlRate, vtmControlModelConfig_.driftLowpassCutoff);

	// Configure intonation.
	eventList_.setMicroIntonation( vtmControlModelConfig_.microIntonation);
	eventList_.setMacroIntonation( vtmControlModelConfig_.macroIntonation);
	eventList_.setSmoothIntonation(vtmControlModelConfig_.smoothIntonation);
	eventList_.setIntonationDrift( vtmControlModelConfig_.intonationDrift);
	eventList_.setRandomIntonation(vtmControlModelConfig_.randomIntonation);
	eventList_.setIntonationFactor(vtmControlModelConfig_.intonationFactor);
}

bool
Controller::nextChunk(const std::string& phoneticString, std::size_t& index, std::size_t& size)
{
	static const std::string token{"/c"};

	const std::size_t startPos = phoneticString.find(token, index);
	if (startPos == std::string::npos) {
		index = phoneticString.size();
		size = 0;
		return false;
	}
	index = startPos;

	const std::size_t endPos = phoneticString.find(token, startPos + token.size());
	if (endPos == std::string::npos) {
		size = phoneticString.size() - startPos;
	} else {
		size = endPos                - startPos;
	}

	for (std::size_t i = index + token.size(), end = index + size; i < end; ++i) {
		if (!std::isspace(phoneticString[i])) {
			LOG_DEBUG("[Controller::nextChunk] Phonetic string chunk: \"" <<
					std::string(phoneticString.begin() + index, phoneticString.begin() + index + size) << '"');
			return true;
		}
	}

	LOG_DEBUG("[Controller::nextChunk] Empty chunk.");

	// The chunk contains only spaces or is empty.
	return false;
}

inline bool ak(std::string s,char c) // अ॒धि॒का॒र॒सू॒च॒कः
{
	return s.find(c)!=std::string::npos;
}
void Controller::vk(std::string s)
{
	{
		std::string ns="";
		for(int i=0;i<s.size();i++)
			if(ak("aAiIuUfFxXeEoOMHkKgGNcCjJYwWqQRtTdDnpPbBmyrlvSzsh",s[i]))ns+=s[i];
		s=ns;
	}
	std::cout<<"\""<<s<<"\""<<std::endl;

	float P[128][22];
	auto set=[this,&P](unsigned char v,std::string name)
	{
		auto PL=model_.postureList().find(name);
		int i;
		for(i=0;i<16;i++)
			P[v][i]=PL->getParameterTarget(i);
		for(i=16;i<22;i++)P[v][i]=0;
		P[v][16]=P[v][12];
	};
	set(' ',"a");
	P[(unsigned char)' '][1]=0;

	set('a', "a");set('A',"ar");set('i', "i");set('I',"ee");set('u',"uu");set('U',"uu");
	set('f',"rr");set('F',"rr");set('x', "l");set('X', "l");set('e', "e");set('E',"er");set('o',"o");set('O',"aw");
	set('M', "m");set('H', "h");
	set('k', "k");set('K', "k");set('g', "g");set('G', "g");set('N',"ng");
	set('c',"ch");set('C',"ch");set('j', "j");set('J', "j");set('Y', "n");
	set('w', "t");set('W', "t");set('q', "d");set('Q', "d");set('R', "t");
	set('t',"th");set('T',"th");set('d',"dh");set('D',"dh");set('n',"th");
	set('p', "p");set('P', "p");set('b', "b");set('B', "b");set('m', "m");
	set('y', "y");set('r',"rr");set('l', "l");set('v', "v");
	set('S', "s");set('z',"sh");set('s', "s");set('h', "h");
	set('V',"f");
	P[(unsigned char)'v'][3]=0;
	P[(unsigned char)'v'][1]=54;
	P[(unsigned char)'v'][14]=0.4;
//	P[(unsigned char)'v'][13]=0.4;

	for(unsigned char i:std::string("tTdDn"))
	{
		P[i][13]=0.1;
		P[i][12]=0.2;
		P[i][3]=0;
	}
	for(unsigned char i:std::string("wWqQR"))
	{
		P[i][16]=0.1;
		P[i][12]=1;
	}
	for(unsigned char i:std::string("Rn"))
	{
		for(int p=0;p<4;p++)
			P[i][p]=P[(unsigned char)'m'][p];
		P[i][15]=P[(unsigned char)'m'][15];
	}
	if(0)for(unsigned char i:std::string("cCjJ"))
	{
		for(int p=3;p<4;p++)
			P[i][p]=P[(unsigned char)'S'][p];
	}
	for(unsigned char i:std::string("KGCJWQTDPB"))
	{
		P[i][3]=14;
	}
	std::vector<float> PL;
	PL.resize(22);
	for(size_t i=0;i<s.size();i++)
	{
		unsigned char v=s[i];
		unsigned char _v=i>0?s[i-1]:' ';
		unsigned char v_=(i<s.size()-1)?s[i+1]:' ';
		if(v=='H')
		{
			if(ak("pP",v_))v='V';
			if(ak("Szs",v_))v=v_;
		}
		float hd=0.15; // ह्र॒स्वा॒व॒धिः॒
		float pd=0.01; // प्र॒थ॒मा॒व॒धिः
		float vd=      // अ॒व॒धिः 
			ak("aiufx",v)?hd
			:ak("AIUFXeEoO",v)?hd*2.0
			:0.1;
		bool ks=ak("aAiIuUfFxXeEoOMgGNjJYqQRdDnbBmyrlv",v); // क॒ण्ठ॒श॒ब्दा॒य
		
		auto ms=[&P](int p,unsigned char v1,unsigned char v2) // म॒ध्य॒स्थि॒तिः
		{
			if(p==2&&ak("KGCJWQTDPB",v1))return (float)40.0;
			else if(ak("aiufxAIUFXeEoO",v1)&&v2=='y')return P[v2][p];
			else if(ak("aiufxAIUFXeEoO",v2)&&v1=='y')return P[v1][p];
			else if(p==15) return std::max(P[v1][p],P[v2][p]);
			else return std::min(P[v1][p],P[v2][p]);
		};
		for(double t=0;t<vd;t+=(float)vtmControlModelConfig_.controlPeriod/1000.0)
		{
			for(int p=0;p<22;p++)
				PL[p]=
					((p>=7&&p<=14)||p==16)?(
					ak("aiufx",v)?
					((t<vd/2.0)?
					(ms(p,_v,v)*(1.0-t*2.0/vd)+P[v][p]*t*2.0/vd)
					:(ms(p,v,v_)*(t*2.0/vd-1)+P[v][p]*(2.0-t*2.0/vd)))
					:ak("AIUFXeEoO",v)?
					((t<hd/2.0)?
					(ms(p,_v,v)*(1.0-t*2.0/hd)+P[v][p]*t*2.0/hd)
					:(vd-t<hd/2.0)?
					(ms(p,v,v_)*(1-(vd-t)*2.0/hd)
					 +P[v][p]*((vd-t)*2.0/hd))
					:P[v][p]
					)					
					:
					(ms(p,_v,v)*(1.0-t/vd)
					 +ms(p,v,v_)*(t/vd))
					)
					:(p==1?
					((t<pd/2.0)?
					(ms(p,_v,v)*(1.0-t*2.0/pd)+P[v][p]*t*2.0/pd)
					:(vd-t<pd/2.0)?
					(ms(p,v,v_)*(1-(vd-t)*2.0/pd)
					 +P[v][p]*((vd-t)*2.0/pd))
					:P[v][p]
					)
					:p==15?
					(ms(p,_v,v)*(1.0-t/vd)
					+ms(p,v,v_)*(t/vd))
					:((t<vd/2.0)?
					(ms(p,_v,v)*(1.0-t*2.0/vd)+P[v][p]*t*2.0/vd)
					:(ms(p,v,v_)*(t*2.0/vd-1)+P[v][p]*(2.0-t*2.0/vd))));
			PL[0]=-10;
			vtmParamList_.push_back(PL);
		}
		
	}
	for(double t=0;t<0.1;t+=(float)vtmControlModelConfig_.controlPeriod/1000.0)
		vtmParamList_.push_back(PL);
	if(false){
	L:
	auto a_PL=model_.postureList().find("a");
	std::vector<float> PL;
	PL.resize(22);
	for(int i=0;i<16;i++)
		PL[i]=a_PL->getParameterTarget(i);
	PL[0]=-10;
	PL[1]=1;
	PL[16]=PL[12];
	PL[17]=0;
	PL[18]=2;
	PL[19]=0;
	PL[20]=0;
	PL[21]=0;

	float ak=0.15,rp=0.05,rk=0.;
	for(int i=0;i<(ak*2+rp*2+rk+0.5)*250;i++)
	{
		PL[17]=0.1+
		((i<ak*250||i>(ak+rp*2+rk)*250)?0
		:(i<(ak+rp)*250)?((float)i-ak*250)/rp/250
		:(i>(ak+rp+rk)*250)?((ak+rp*2.0+rk)*250-i)/rp/250
		:1)/2.0;

		if(i<ak*250*0.5){PL[1]=60.0*(float)i/ak/0.5/250;}
		else if((ak*2+rp*2+rk)*250-i<ak*0.25*250)PL[1]=60.0*std::max(((ak*2+rp*2+rk)*250-i)/ak/0.25/250,0.0);
		else PL[1]=60-PL[17]*2*5;

		PL[16]=1.8-PL[17];
		PL[18]=PL[17];
		//PL[11]=1.8-PL[17];
		PL[12]=a_PL->getParameterTarget(12)+(0.2-a_PL->getParameterTarget(12))*PL[17]*2.0;
		vtmParamList_.push_back(PL);
	}
	}
}
void
Controller::getParametersFromPhoneticString(const std::string& phoneticString)
{
	vtmParamList_.clear();
	initUtterance();

	if(sk)vk(phoneticString);
	else if (vtmControlModelConfig_.phoStrFormat == PhoneticStringFormat::mbrola) {
		if (!pho1Parser_) {
			pho1Parser_ = std::make_unique<Pho1Parser>(configDirPath_.c_str(), model_, eventList_);
		}

		eventList_.setMicroIntonation(true);
		eventList_.setMacroIntonation(true);
		eventList_.setSmoothIntonation(false);

		eventList_.setUp();
		pho1Parser_->parse(phoneticString);
		eventList_.generateOutput(vtmParamList_);
	} else {
		if (!phoneticStringParser_) {
			phoneticStringParser_ = std::make_unique<PhoneticStringParser>(configDirPath_.c_str(), model_, eventList_);
		}

		std::size_t index = 0, size = 0;
		while (index < phoneticString.size()) {
			if (nextChunk(phoneticString, index, size)) {
				eventList_.setUp();

				phoneticStringParser_->parse(&phoneticString[index], size);

				eventList_.generateEventList();
				eventList_.applyIntonation();
				eventList_.generateOutput(vtmParamList_);
			}

			index += size;
		}
	}
}

void
Controller::getParametersFromEventList()
{
	vtmParamList_.clear();

	initUtterance();

	eventList_.clearMacroIntonation();
	eventList_.prepareMacroIntonationInterpolation();
	eventList_.generateOutput(vtmParamList_);
}

void
Controller::getParametersFromStream(std::istream& in)
{
	vtmParamList_.clear();
	std::string line;
	const std::size_t numParam = model_.parameterList().size();
	std::vector<float> param(numParam);

	unsigned int lineNumber = 1;
	while (std::getline(in, line)) {
		std::istringstream lineStream(line);

		for (std::size_t i = 0; i < numParam; ++i) {
			lineStream >> param[i];
		}
		if (!lineStream) {
			THROW_EXCEPTION(VTMException, "Could not read vocal tract parameters from stream (line number " << lineNumber << ").");
		}

		vtmParamList_.push_back(param);
		++lineNumber;
	}
}

void
Controller::synthesizePhoneticStringToFile(const std::string& phoneticString, const char* vtmParamFile, const char* outputFile)				
{
	getParametersFromPhoneticString(phoneticString);
	if (vtmParamFile) writeVTMParameterFile(vtmParamList_, vtmParamFile);
	synthesizeToFile(outputFile);
}

void
Controller::synthesizePhoneticStringToBuffer(const std::string& phoneticString, const char* vtmParamFile, std::vector<float>& buffer)
{
	getParametersFromPhoneticString(phoneticString);
	if (vtmParamFile) writeVTMParameterFile(vtmParamList_, vtmParamFile);
	synthesizeToBuffer(buffer);
}

void
Controller::synthesizeFromEventListToFile(const char* vtmParamFile, const char* outputFile)
{
	getParametersFromEventList();
	if (vtmParamFile) writeVTMParameterFile(vtmParamList_, vtmParamFile);
	synthesizeToFile(outputFile);
}

void
Controller::synthesizeFromEventListToBuffer(const char* vtmParamFile, std::vector<float>& buffer)
{
	getParametersFromEventList();
	if (vtmParamFile) writeVTMParameterFile(vtmParamList_, vtmParamFile);
	synthesizeToBuffer(buffer);
}

void
Controller::synthesizeToFile(const char* outputFile)
{
	if (!outputFile) return;

	if (!vtm_->outputBuffer().empty()) vtm_->reset();
	synthesize(vtmParamList_);
	vtm_->finishSynthesis();
	writeOutputToFile(outputFile, outputScale_);
}

void
Controller::synthesizeToFile(std::vector<std::vector<float>>& vtmParamList, const char* vtmParamFile, const char* outputFile)
{
	if (vtmParamFile) writeVTMParameterFile(vtmParamList, vtmParamFile);

	if (!vtm_->outputBuffer().empty()) vtm_->reset();
	synthesize(vtmParamList);
	vtm_->finishSynthesis();
	float scale;
	writeOutputToFile(outputFile, scale);
}

void
Controller::synthesizeToBuffer(std::vector<float>& outputBuffer)
{
	if (!vtm_->outputBuffer().empty()) vtm_->reset();
	synthesize(vtmParamList_);
	vtm_->finishSynthesis();
	writeOutputToBuffer(outputBuffer, outputScale_);
}

void
Controller::synthesizeToBuffer(std::vector<std::vector<float>>& vtmParamList, const char* vtmParamFile, std::vector<float>& outputBuffer)
{
	if (vtmParamFile) writeVTMParameterFile(vtmParamList, vtmParamFile);

	if (!vtm_->outputBuffer().empty()) vtm_->reset();
	synthesize(vtmParamList);
	vtm_->finishSynthesis();
	float scale;
	writeOutputToBuffer(outputBuffer, scale);
}

void
Controller::synthesizeToFile(std::istream& inputStream, const char* outputFile)
{
	getParametersFromStream(inputStream);
	synthesizeToFile(outputFile);
}

void
Controller::synthesize(std::vector<std::vector<float>>& vtmParamList)
{
	if (vtmParamList.empty()) return;

	// Duplicate the last set of parameters, to help the interpolation.
	vtmParamList.push_back(vtmParamList.back());

	// Number of internal sample rate periods in each control rate period.
	const unsigned int controlSteps = static_cast<unsigned int>(std::rint(vtm_->internalSampleRate() / vtmControlModelConfig_.controlRate));
	const float coef = 1.0f / controlSteps;

	const std::size_t numParam = sk?22:model_.parameterList().size();
	std::vector<float> currentParameter(numParam);
	std::vector<float> currentParameterDelta(numParam);

	// For each control period:
	for (std::size_t i = 1, size = vtmParamList.size(); i < size; ++i) {
		// Calculates the current parameter values, and their
		// associated sample-to-sample delta values.
		for (std::size_t j = 0; j < numParam; ++j) {
			currentParameter[j] = vtmParamList[i - 1][j];
			currentParameterDelta[j] = (vtmParamList[i][j] - currentParameter[j]) * coef;
		}

		// For each step in a control period:
		for (std::size_t j = 0; j < controlSteps; ++j) {
			vtm_->setAllParameters(currentParameter);
			vtm_->execSynthesisStep();

			// Do linear interpolation.
			for (std::size_t k = 0; k < numParam; ++k) {
				currentParameter[k] += currentParameterDelta[k];
			}
		}
	}
}

void
Controller::writeOutputToFile(const char* outputFile, float& scale)
{
	if (!outputFile) {
		THROW_EXCEPTION(MissingValueException, "Missing output file name.");
	}
	const std::vector<float>& audioData = vtm_->outputBuffer();
	WAVEFileWriter fileWriter(outputFile, 1, audioData.size(), vtm_->outputSampleRate());

	scale = VTM::Util::calculateOutputScale(audioData);
	for (std::size_t i = 0, end = audioData.size(); i < end; ++i) {
		fileWriter.writeSample(audioData[i] * scale);
	}
}

void
Controller::writeOutputToBuffer(std::vector<float>& outputBuffer, float& scale)
{
	const std::vector<float>& audioData = vtm_->outputBuffer();
	outputBuffer.resize(audioData.size());

	scale = VTM::Util::calculateOutputScale(audioData);
	for (std::size_t i = 0, end = audioData.size(); i < end; ++i) {
		outputBuffer[i] = audioData[i] * scale;
	}
}

void
Controller::writeVTMParameterFile(const std::vector<std::vector<float>>& vtmParamList, const char* vtmParamFile)
{
	if (!vtmParamFile) {
		THROW_EXCEPTION(MissingValueException, "Missing output VTM parameter file name.");
	}
	std::ofstream out(vtmParamFile, std::ios_base::binary);
	if (!out) {
		THROW_EXCEPTION(IOException, "Could not open the file " << vtmParamFile << '.');
	}

	for (auto& param : vtmParamList) {
		if (param.empty()) {
			THROW_EXCEPTION(InvalidValueException, "Empty parameter set.");
		}
		out << param[0];
		for (unsigned int i = 1, size = param.size(); i < size; ++i) {
			out << ' ' << param[i];
		}
		out << '\n';
	}
}

} /* namespace VTMControlModel */
} /* namespace GS */
