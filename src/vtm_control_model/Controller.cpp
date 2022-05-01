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
	if(0){
	auto a_PL=model_.postureList().find("i");
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
	{
		std::string ns="";
		for(size_t i=0;i<s.size();i++)
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
		{
			P[v][i]=PL->getParameterTarget(i);
			if(false&&i>=7&&i<=14&&P[v][i]<0.11)P[v][i]=0;
		}
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
	set('y', "y");set('r', "t");set('l', "l");set('v', "u");
	set('S', "s");set('z',"sh");set('s', "s");set('h', "h");
	set('V',"f");
	P[(unsigned char)'v'][3]=0;
	P[(unsigned char)'v'][1]=54;
	P[(unsigned char)'v'][14]=0.4;
//	P[(unsigned char)'v'][13]=0.4;

	for(unsigned char i:std::string("tTdDn"))
	{
		P[i][13]=0.1;
		P[i][12]=0.3;
	//	P[i][16]=0.2;
	//	P[i][11]=0.3;
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
	for(unsigned char i:std::string("jJ"))
	{
		for(int p=3;p<4;p++)
			P[i][p]=40;
	}
	for(unsigned char i:std::string("cC"))
	{
		for(int p=3;p<4;p++)
			P[i][p]=18;
	}
	for(unsigned char i:std::string("KGCJWQTDPB"))
	{
		P[i][2]=14;
	}
	for(int i=4;i<=6;i++)
		P[(unsigned char)'S'][i]=P[(unsigned char)'c'][i];
	for(int i=7;i<=16;i++)
	{
		float nt=10;
		for(unsigned char j:std::string("aium"))
			nt=std::min(nt,P[j][i]);
		P[(unsigned char)'M'][i]=nt;
	}
	for(int i=16;i<=18;i++)
	{
		//if(i==17||i==18)
		//	P[(unsigned char)'r'][i]=0.6;
		//else
		if(i==16)P[(unsigned char)'r'][i]=1.8;
		if(i==17||i==18)P[(unsigned char)'r'][i]=0.1;
		else P[(unsigned char)'r'][i]=2;
	}
	P[(unsigned char)'r'][12]=0.47;
	P[(unsigned char)'r'][1]=55;
	P[(unsigned char)'M'][15]=1.0;
	P[(unsigned char)'z'][6]=1700;
	std::vector<float> PL;
	PL.resize(22);
	for(size_t i=0;i<s.size();i++)
	{
		unsigned char v=s[i];
		unsigned char pv=i>0?s[i-1]:' ';
		unsigned char vc=(i<s.size()-1)?s[i+1]:' ';
		if(v=='H')
		{
			if(ak("pP",vc))v='V';
			if(ak("Szs",vc))v=vc;
		}
		float vnd=0.1*0.75;
		float hd=0.15*0.75; // ह्र॒स्वा॒व॒धिः॒
		float pd=-0.01; // प्र॒थ॒मा॒व॒धिः
		auto vdd=[hd,vnd,pd](unsigned char v)
		{
			return ak("aiufx",v)?hd
			:ak("AIUFXeEoO",v)?hd*2.0
			:vnd;
		};
		float vd=vdd(v);     // अ॒व॒धिः 
		if(ak("aiufx",pv))
		{	
			float mk=hd+vd;
			int vs=1;
			for(size_t j=i+1;j<s.size();j++)
			{
				if(ak("kKgGNcCjJYwWqQRtTdDnpPbBmyrlvSzsh",s[j]))
				{
					mk+=vdd(s[j]);
					vs++;
				}
				else break;
			}
			if(vs>1)vd+=std::max(3*hd-mk,(float)0.0);
		}
		[[maybe_unused]] bool ks=ak("aAiIuUfFxXeEoOMgGNjJYqQRdDnbBmyrlv",v); // क॒ण्ठ॒श॒ब्दाय॑
		
		auto ms=[&P](int p,unsigned char v1,unsigned char v2) // म॒ध्य॒स्थि॒तिः
		{
			float ykms=0.5;
			if(ak(" ",v1))return P[v2][p];
			else if(ak(" ",v2))return P[v1][p];
			else if(p==2&&ak("KGCJWQTDPB",v1))return (float)60.0; // unused
			else if(ak("aiufxAIUFXeEoO",v1)&&v2=='y')return (P[v1][p]*((float)1.0-ykms)+P[v2][p]*ykms);
			else if(ak("aiufxAIUFXeEoO",v2)&&v1=='y')return (P[v1][p]*ykms+P[v2][p]*((float)1.0-ykms));
			else if(p==15) return std::max(P[v1][p],P[v2][p]);
			else if(p==17||p==18)return std::max(P[v1][p],P[v2][p]);
			else if((v1=='r'||v2=='r')&&p==16&&v1!='R'&&v2!='R')return P[(unsigned char)'r'][p];
			else return std::min(P[v1][p],P[v2][p]);
		};
		auto ayauk=[&P](int p,unsigned char v,int k)
		{
			if(v=='E')return k?P[(unsigned char)'i'][p]:P[(unsigned char)'a'][p];
			else if(v=='O')return k?P[(unsigned char)'u'][p]:P[(unsigned char)'a'][p];
			else return P[v][p];
		};
		//TODO: र॒का॒रः।  अ॒नु॒स्वा॒रः।  ह्र॒स्वा॒र्द्ध॒का॒ले आ॑स्यपरि॒वर्त्त॑नम्।
		for(double t=0;t<vd;t+=(float)vtmControlModelConfig_.controlPeriod/1000.0)
		{
			float mpk=hd/2.0; // म॒हा॒प्रा॒ण॒का॒लः
			int p=0;
			PL[p]=-6.0*(1.0+0.05*(float)rand()/(float)RAND_MAX);
			{
				p=1;
				float nd=std::min(0.03,hd*0.2);
				if(ak("KGCJWQTDPB",pv)&&t<mpk)PL[p]=P[v][p]*t/mpk;
				else if((ak(" ",vc)||(ak("aiufxAIUFXeEoO",vc)&&ak("aiufxAIUFXeEoO",v)))&&vd-t<nd)
					PL[p]=P[v][p]*((vd-t)/nd);
				else if((ak(" ",pv)||(ak("aiufxAIUFXeEoO",pv)&&ak("aiufxAIUFXeEoO",v)))&&t<nd)
						PL[p]=P[v][p]*t/nd;
				else PL[p]=(t<pd/2.0)?
						(ms(p,pv,v)*(1.0-t*2.0/pd)+P[v][p]*t*2.0/pd)
						:(vd-t<pd/2.0)?
						(ms(p,v,vc)*(1-(vd-t)*2.0/pd)
						 +P[v][p]*((vd-t)*2.0/pd))
						:P[v][p];
			}
			p=2;
			if(ak("KGCJWQTDPB",pv)&&t<mpk)PL[p]=P[pv][p];
			else PL[p]=P[v][p];
			p=3;
			if(ak("cCjJ",v))
			{
				PL[p]=ak("cCjJ",vc)?0:t<vd/2.0?0:P[v][p];	
			}
			else PL[p]=P[v][p];
			for(p=4;p<7;p++)
				PL[p]=P[v][p];
			for(p=7;p<22;p++)
			{
				if(p!=15)PL[p]=
					(ak("aiufx",v)||ak("kKgGNcCjJYwWqQRtTdDnpPbBmyrlvSzsh",v))?
					((t<vd/2.0)?
					(ms(p,pv,v)*(1.0-t*2.0/vd)+P[v][p]*t*2.0/vd)
					:(ms(p,v,vc)*(t*2.0/vd-1)+P[v][p]*(2.0-t*2.0/vd)))
					:ak("AIUFXeEoO",v)?
					((t<hd/2.0)?
					(ms(p,pv,v)*(1.0-t*2.0/hd)
					 +(ak("EO",v)?ayauk(p,v,0)
					 :(v=='A'&&pv=='n')?P[(unsigned char)'a'][p]// ना॒ना॒कृ॒ते।
					 :P[v][p])*t*2.0/hd)
					:(vd-t<hd/2.0)?
					(ms(p,v,vc)*(1-(vd-t)*2.0/hd)
					 +(ak("EO",v)?ayauk(p,v,1):P[v][p])*((vd-t)*2.0/hd))
					:(ak("EO",v)?(ayauk(p,v,0)*(vd-hd/2.0-t)/(vd-hd)+ayauk(p,v,1)*(t-hd/2.0)/(vd-hd))
						// ना॒ना॒कृ॒ते।
					 :(v=='A'&&pv=='n')?P[(unsigned char)'a'][p]*(vd-hd/2.0-t)/(vd-hd)+P[v][p]*(t-hd/2.0)/(vd-hd)
					 :P[v][p])
					)
					:
					(ms(p,pv,v)*(1.0-t/vd)
					 +ms(p,v,vc)*(t/vd))
					;
			}
			if(v=='r')
			{
				PL[17]=P[v][17]+(t<vd/2.0?t*2.0/vd:(2.0-t*2.0/vd))/2.0;
				PL[18]=PL[17];
				PL[16]=P[v][16]-PL[17];
				PL[12]=((t<vd/2.0)?ms(12,pv,v):ms(12,v,vc))*(1-PL[17]*2.0)+0.2*PL[17]*2.0;
			}
			{
				p=15;
				double nd=std::min(hd,vd)/2.0;
				PL[p]=true?(ms(p,pv,v)*(1.0-t/vd)
					+ms(p,v,vc)*(t/vd)):
					((t<nd)?
					(ms(p,pv,v)*(1.0-t/nd)+P[v][p]*t/nd)
					:(vd-t<nd)?
					(ms(p,v,vc)*(1-(vd-t)/nd)+P[v][p]*((vd-t)/nd))
					:P[v][p]);
			}
		
			vtmParamList_.push_back(PL);
		}
		if(ak("aiufxAIUFXeEoO",vc)&&ak("aiufxAIUFXeEoO",v))
			for(double t=0;t<hd*0.3;t+=(float)vtmControlModelConfig_.controlPeriod/1000.0)
				vtmParamList_.push_back(PL);

	}
	for(double t=0;t<0.1;t+=(float)vtmControlModelConfig_.controlPeriod/1000.0)
		vtmParamList_.push_back(PL);
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
