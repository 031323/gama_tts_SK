/***************************************************************************
 *  Copyright 2014 Marcelo Y. Matuda                                       *
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
// This file was created by Marcelo Y. Matuda, and code/information
// from Gnuspeech was added to it later.

#ifndef TRM_CONTROL_MODEL_TRANSITION_H_
#define TRM_CONTROL_MODEL_TRANSITION_H_

#include <string>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Equation.h"
#include "Exception.h"



namespace GS {
namespace TRMControlModel {

class Model;

class Transition {
public:
	enum Type {
		TYPE_INVALID    = 0,
		TYPE_DIPHONE    = 2,
		TYPE_TRIPHONE   = 3,
		TYPE_TETRAPHONE = 4
	};

	struct PointOrSlope {
		virtual ~PointOrSlope() {}
		virtual bool isSlopeRatio() const = 0;
	};
	typedef std::unique_ptr<PointOrSlope> PointOrSlope_ptr;

	struct Point : PointOrSlope {
		enum Type {
			TYPE_INVALID    = 0,
			TYPE_DIPHONE    = 2, // Posture 1 --> Posture 2
			TYPE_TRIPHONE   = 3, // Posture 2 --> Posture 3
			TYPE_TETRAPHONE = 4  // Posture 3 --> Posture 4
		};
		Type type;
		float value;
		bool isPhantom;

		// If timeExpression is not empty, time = timeExpression, otherwise time = freeTime.
		std::string timeExpression;
		float freeTime; // milliseconds

		Point() : type(TYPE_INVALID), value(0.0), isPhantom(false), timeExpression(), freeTime(0.0) {}
		virtual ~Point() {}

		virtual bool isSlopeRatio() const { return false; }

		static Type getTypeFromName(const std::string& typeName) {
			if (typeName == "diphone") {
				return Transition::Point::TYPE_DIPHONE;
			} else if (typeName == "triphone") {
				return Transition::Point::TYPE_TRIPHONE;
			} else if (typeName == "tetraphone") {
				return Transition::Point::TYPE_TETRAPHONE;
			} else {
				THROW_EXCEPTION(TRMControlModelException, "Invalid transition point type: " << typeName << '.');
			}
		}
	private:
		Point(const Point&);
		Point& operator=(const Point&);
	};
	typedef std::unique_ptr<Point> Point_ptr;

	struct Slope {
		float slope;
		float displayTime;

		Slope() : slope(0.0), displayTime(0.0) {}
	private:
		Slope(const Slope&);
		Slope& operator=(const Slope&);
	};
	typedef std::unique_ptr<Slope> Slope_ptr;

	struct SlopeRatio : PointOrSlope {
		std::vector<Point_ptr> pointList;
		std::vector<Slope_ptr> slopeList;

		SlopeRatio() {}
		virtual ~SlopeRatio() {}

		virtual bool isSlopeRatio() const { return true; }
		double totalSlopeUnits() const;
	private:
		SlopeRatio(const SlopeRatio&);
		SlopeRatio& operator=(const SlopeRatio&);
	};

	Transition(
		const std::string& groupName,
		const std::string& name,
		Type type,
		bool special)
			: groupName_(groupName)
			, name_(name)
			, type_(type)
			, special_(special)
	{
	}

	const std::string& name() const {
		return name_;
	}
	Type type() const {
		return type_;
	}
	const std::string& groupName() const {
		return groupName_;
	}
	bool special() const {
		return special_;
	}

	std::vector<PointOrSlope_ptr>& pointOrSlopeList() { return pointOrSlopeList_; }
	const std::vector<PointOrSlope_ptr>& pointOrSlopeList() const { return pointOrSlopeList_; }

	static double getPointTime(const Transition::Point& point, const Model& model);
	static void getPointData(const Transition::Point& point, const Model& model,
					double& time, double& value);
	static void getPointData(const Transition::Point& point, const Model& model,
					double baseline, double delta, double min, double max,
					double& time, double& value);
	static Type getTypeFromName(const std::string& typeName) {
		if (typeName == "diphone") {
			return Transition::TYPE_DIPHONE;
		} else if (typeName == "triphone") {
			return Transition::TYPE_TRIPHONE;
		} else if (typeName == "tetraphone") {
			return Transition::TYPE_TETRAPHONE;
		} else {
			THROW_EXCEPTION(TRMControlModelException, "Invalid transition type: " << typeName << '.');
		}
	}
private:
	std::string groupName_;
	std::string name_;
	Type type_;
	bool special_;
	std::vector<PointOrSlope_ptr> pointOrSlopeList_;
};

} /* namespace TRMControlModel */
} /* namespace GS */

#endif /* TRM_CONTROL_MODEL_TRANSITION_H_ */
