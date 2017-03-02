#include "Curve.h"
#include <fstream>
#include <iostream>
#include <cmath>

const std::string gTypeKey = "Type";
const std::string gAnchorsKey = "Anchors";
const std::string gAnchorPosKey = "Pos";
const std::string gSegmentDurationKey = "SegmentDuration";

cCurve::tAnchor::tAnchor()
{
	mPos.setZero();
	mTangent.setZero();
}

cCurve::tAnchor::~tAnchor()
{
}

cCurve::eCurveType cCurve::ParseCurveType(const std::string& str)
{
	// convert a string into the corresponding curve type
	eCurveType curve_type = eCurveTypeCatmullRom;

	if (str == "catmull_rom")
	{
		curve_type = eCurveTypeCatmullRom;
	}
	else if (str == "b_spline")
	{
		curve_type = eCurveTypeBSpline;
	}
	else
	{
		printf("Unsupported curve type %s\n", str.c_str());
		assert(false); // unsupported curve type
	}

	return curve_type;
}

cCurve::cCurve()
{
	mSegmentDuration = 1;
	mCurveType = eCurveTypeCatmullRom;
}

cCurve::~cCurve()
{
}

bool cCurve::Load(const std::string& file)
{
	// load a set of anchor points and other parameters from a text file
	// file should be formatted as a JSON
	bool succ = true;
	Clear();

	std::ifstream f_stream(file.c_str());
	Json::Value root;
	Json::Reader reader;
	succ = reader.parse(f_stream, root);
	f_stream.close();

	if (succ)
	{
		// parse misc parameters from the file
		std::string type_str = root.get(gTypeKey, "").asString();
		mCurveType = ParseCurveType(type_str);
		mSegmentDuration = root.get(gSegmentDurationKey, 1).asDouble();

		// parse the list of anchors
		if (!root[gAnchorsKey].isNull())
		{
			auto anchors_json = root.get(gAnchorsKey, 0);
			succ &= ParseAnchors(anchors_json);
		}
	}

	if (succ)
	{
		// provides some examples of how to work with the anchor data structure
		PrintAnchors();
	}
	else
	{
		Clear();
	}

	return succ;
}

void cCurve::Clear()
{
	mAnchors.clear();
	mSegmentDuration = 1;
}

int cCurve::GetNumAnchors() const
{
	return static_cast<int>(mAnchors.size());
}

const Eigen::VectorXd& cCurve::GetAnchorPos(int i) const
{
	return mAnchors[i].mPos;
}

const Eigen::VectorXd& cCurve::GetAnchorTangent(int i) const
{
	return mAnchors[i].mTangent;
}

int cCurve::GetNumSegments() const
{
	// computes the number of curve segments from the number of anchors
	int num_anchors = GetNumAnchors();
	int num_segs = 0;
	switch (mCurveType)
	{
	case eCurveTypeCatmullRom:
		num_segs = num_anchors - 1;
		break;
	case eCurveTypeBSpline:
		num_segs = num_anchors + 1;
		break;
	default:
		assert(false); // unsuppoted curve type
		break;
	}

	return num_segs;
}

int cCurve::GetDim() const
{
	// dimension of each anchor
	int dim = 0;
	if (GetNumAnchors() > 0)

	{
		const auto& pos = GetAnchorPos(0);
		dim = pos.size();
	}
	return dim;
}

void cCurve::Eval(double time, Eigen::VectorXd& out_result) const
{
	// TODO (CPSC426): Evaluates a parametric curve at the given time
	out_result = Eigen::VectorXd::Zero(GetDim());
	//out_result[0] = time; // stub, position just moves along the x-axis with time
	

	// first build the basis matrix M for the current curve type (mCurveType)

	Eigen::Matrix4d M_matrix;
	double coefficient;
	switch (mCurveType) { // set M matrix and coefficient (1/2 or 1/6) depending on the curve type
	case eCurveTypeCatmullRom:
		coefficient = 0.5;
		M_matrix << -1.0, 3.0, -3.0, 1.0,
			2.0, -5.0, 4.0, -1.0,
			-1.0, 0.0, 1.0, 0.0,
			0.0, 2.0, 0.0, 0.0;
		break;
	case eCurveTypeBSpline:
		coefficient = 1.0 / 6.0;
		M_matrix << -1.0, 3.0, -3.0, 1.0,
			3.0, -6.0, 3.0, 0.0,
			-3.0, 0.0, 3.0, 0.0,
			1.0, 4.0, 1.0, 0.0;
		break;
	default:
		assert(false);
		break;
	}

	M_matrix = coefficient * M_matrix;
	// then build the T polynomial vector
	int segment = (int)time;
	double t = time - segment; // this is how we clamp time

	Eigen::Vector4d T_vector;
	T_vector << pow(t, 3), pow(t, 2), t, 1.0;

	// finally build the geometry matrix G for the curve segment
	int dim = GetDim();
	int anc = GetNumAnchors();
	int seg = GetNumSegments();
	Eigen::MatrixXd G_matrix(anc, dim);

	for (int m = 0; m < anc; ++m) {
		const tAnchor &anchor = mAnchors[m];
		for (int n = 0; n < dim; ++n) {
			G_matrix(m, n) = anchor.mPos[n];
		}
	}

	Eigen::MatrixXd beg(1, dim);
	Eigen::MatrixXd end(1, dim);
	beg << G_matrix.row(0);
	end << G_matrix.row(anc - 1);

	int anchor_beg, anchor_end;
	GetAnchors(segment, anchor_beg, anchor_end);

	Eigen::MatrixXd G_matrixCatmull_beg(anc + 1, dim);
	Eigen::MatrixXd G_matrixCatmull_end(anc + 1, dim);
	Eigen::MatrixXd G_matrixBSpline_beg(anc + 4, dim);
	Eigen::MatrixXd G_matrixBSpline_end(anc + 2, dim);

	switch (mCurveType) {
	case eCurveTypeCatmullRom:
	
		G_matrixCatmull_beg << beg, G_matrix;
		G_matrixCatmull_end << G_matrix, end;

		if (segment == 0) {
			out_result = T_vector.transpose() * M_matrix * G_matrixCatmull_beg.block(0, 0, 4, dim);
		}
		else if ((0 < segment) && (segment < seg - 1)) {
			out_result = T_vector.transpose() * M_matrix * G_matrix.block(std::max(anchor_beg, 0), 0, 4, dim);
		}
		else if (seg == seg - 1) {
			out_result = T_vector.transpose() * M_matrix * G_matrixCatmull_end.block(anchor_beg, 0, 4, dim);
		}
		break;
	case eCurveTypeBSpline:
		
		G_matrixBSpline_beg << beg, beg, G_matrix, end, end;
		G_matrixBSpline_end << G_matrix, end;
		if (segment == 0) {
			out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_beg.block(0, 0, 4, dim);
		}
		else if ((0 < segment) && (segment < seg - 1)) {
			if (anchor_beg < 0) {
				out_result = T_vector.transpose() * M_matrix *G_matrixBSpline_beg.block(anchor_beg + 2, 0, 4, dim);
			}
			else {
				out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_beg.block(anchor_beg, 0, 4, dim);
			}
		}
		else if (segment == seg - 1) {
			out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_end.block(anchor_beg, 0, 4, dim);
		}
		break;
	}

}

void cCurve::EvalTangent(double time, Eigen::VectorXd& out_result) const
{
	// TODO (CPSC426): Evaluates the first derivative of a curve
	out_result = Eigen::VectorXd::Zero(GetDim()); // stub
	

	// so this is going to be the same as Eval, except the time is different. we're using the first derivative of t.

	// first build the basis matrix M for the current curve type (mCurveType)

	Eigen::Matrix4d M_matrix;
	double coefficient;
	switch (mCurveType) { // set M matrix and coefficient (1/2 or 1/6) depending on the curve type
	case eCurveTypeCatmullRom:
		coefficient = 0.5;
		M_matrix << -1.0, 3.0, -3.0, 1.0,
			2.0, -5.0, 4.0, -1.0,
			-1.0, 0.0, 1.0, 0.0,
			0.0, 2.0, 0.0, 0.0;
		break;
	case eCurveTypeBSpline:
		coefficient = 1.0 / 6.0;
		M_matrix << -1.0, 3.0, -3.0, 1.0,
			3.0, -6.0, 3.0, 0.0,
			-3.0, 0.0, 3.0, 0.0,
			1.0, 4.0, 1.0, 0.0;
		break;
	default:
		assert(false);
		break;
	}
	M_matrix = coefficient * M_matrix;
	// then build the T polynomial vector
	int segment = (int)time;
	double t = time - segment;
	Eigen::Vector4d T_vector;
	T_vector << 3.0 *pow(t, 2), 2.0 * t, 1.0, 0.0; // in this case, we use the first derivatives... pow(t,3) -> 3*pow(t,2), etc.

	// finally build the geometry matrix G for the curve segment
	int dim = GetDim();
	int anc = GetNumAnchors();
	int seg = GetNumSegments();
	Eigen::MatrixXd G_matrix(anc, dim);

	for (int m = 0; m < anc; ++m) {
		const tAnchor &anchor = mAnchors[m];
		for (int n = 0; n < dim; ++n) {
			G_matrix(m, n) = anchor.mPos[n];
		}
	}

	Eigen::MatrixXd beg(1, dim);
	Eigen::MatrixXd end(1, dim);
	beg << G_matrix.row(0);
	end << G_matrix.row(anc - 1);

	int anchor_beg, anchor_end;
	GetAnchors(segment, anchor_beg, anchor_end);

	Eigen::MatrixXd G_matrixCatmull_beg(anc + 1, dim);
	Eigen::MatrixXd G_matrixCatmull_end(anc + 1, dim);
	Eigen::MatrixXd G_matrixBSpline_beg(anc + 4, dim);
	Eigen::MatrixXd G_matrixBSpline_end(anc + 2, dim);

	switch (mCurveType) {
	case eCurveTypeCatmullRom:

		G_matrixCatmull_beg << beg, G_matrix;
		G_matrixCatmull_end << G_matrix, end;

		if (segment == 0) {
			out_result = T_vector.transpose() * M_matrix * G_matrixCatmull_beg.block(0, 0, 4, dim);
		}
		else if ((0 < segment) && (segment < seg - 1)) {
			out_result = T_vector.transpose() * M_matrix * G_matrix.block(std::max(anchor_beg, 0), 0, 4, dim);
		}
		else if (seg == seg - 1) {
			out_result = T_vector.transpose() * M_matrix * G_matrixCatmull_end.block(anchor_beg, 0, 4, dim);
		}
		break;
	case eCurveTypeBSpline:

		G_matrixBSpline_beg << beg, beg, G_matrix, end, end;
		G_matrixBSpline_end << G_matrix, end;
		if (segment == 0) {
			out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_beg.block(0, 0, 4, dim);
		}
		else if ((0 < segment) && (segment < seg - 1)) {
			if (anchor_beg < 0) {
				out_result = T_vector.transpose() * M_matrix *G_matrixBSpline_beg.block(anchor_beg + 2, 0, 4, dim);
			}
			else {
				out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_beg.block(anchor_beg, 0, 4, dim);
			}
		}
		else if (segment == seg - 1) {
			out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_end.block(anchor_beg, 0, 4, dim);
		}
		break;
	}
}

void cCurve::EvalNormal(double time, Eigen::VectorXd& out_result) const
{
	// TODO (CPSC426): Evaluates the second derivative of a curve
	out_result = Eigen::VectorXd::Zero(GetDim()); // stub

	//again, this is the same!! second derivatives!!


	// first build the basis matrix M for the current curve type (mCurveType)

	Eigen::Matrix4d M_matrix;
	double coefficient;
	switch (mCurveType) { // set M matrix and coefficient (1/2 or 1/6) depending on the curve type
	case eCurveTypeCatmullRom:
		coefficient = 0.5;
		M_matrix << -1.0, 3.0, -3.0, 1.0,
			2.0, -5.0, 4.0, -1.0,
			-1.0, 0.0, 1.0, 0.0,
			0.0, 2.0, 0.0, 0.0;
		break;
	case eCurveTypeBSpline:
		coefficient = 1.0 / 6.0;
		M_matrix << -1.0, 3.0, -3.0, 1.0,
			3.0, -6.0, 3.0, 0.0,
			-3.0, 0.0, 3.0, 0.0,
			1.0, 4.0, 1.0, 0.0;
		break;
	default:
		assert(false);
		break;
	}

	M_matrix = coefficient * M_matrix;
	// then build the T polynomial vector
	int segment = (int)time;
	double t = time - segment;
	Eigen::Vector4d T_vector;
	T_vector << 6.0 * t, 2.0, 0.0, 0.0; // in this case, we use the second derivatives... pow(t,3) -> 3*pow(t,2) -> 6*t, etc.

	// finally build the geometry matrix G for the curve segment
	int dim = GetDim();
	int anc = GetNumAnchors();
	int seg = GetNumSegments();
	Eigen::MatrixXd G_matrix(anc, dim);

	for (int m = 0; m < anc; ++m) {
		const tAnchor &anchor = mAnchors[m];
		for (int n = 0; n < dim; ++n) {
			G_matrix(m, n) = anchor.mPos[n];
		}
	}

	Eigen::MatrixXd beg(1, dim);
	Eigen::MatrixXd end(1, dim);
	beg << G_matrix.row(0);
	end << G_matrix.row(anc - 1);

	int anchor_beg, anchor_end;
	GetAnchors(segment, anchor_beg, anchor_end);

	Eigen::MatrixXd G_matrixCatmull_beg(anc + 1, dim);
	Eigen::MatrixXd G_matrixCatmull_end(anc + 1, dim);
	Eigen::MatrixXd G_matrixBSpline_beg(anc + 4, dim);
	Eigen::MatrixXd G_matrixBSpline_end(anc + 2, dim);

	switch (mCurveType) {
	case eCurveTypeCatmullRom:

		G_matrixCatmull_beg << beg, G_matrix;
		G_matrixCatmull_end << G_matrix, end;

		if (segment == 0) {
			out_result = T_vector.transpose() * M_matrix * G_matrixCatmull_beg.block(0, 0, 4, dim);
		}
		else if ((0 < segment) && (segment < seg - 1)) {
			out_result = T_vector.transpose() * M_matrix * G_matrix.block(std::max(anchor_beg, 0), 0, 4, dim);
		}
		else if (seg == seg - 1) {
			out_result = T_vector.transpose() * M_matrix * G_matrixCatmull_end.block(anchor_beg, 0, 4, dim);
		}
		break;
	case eCurveTypeBSpline:

		G_matrixBSpline_beg << beg, beg, G_matrix, end, end;
		G_matrixBSpline_end << G_matrix, end;
		if (segment == 0) {
			out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_beg.block(0, 0, 4, dim);
		}
		else if ((0 < segment) && (segment < seg - 1)) {
			if (anchor_beg < 0) {
				out_result = T_vector.transpose() * M_matrix *G_matrixBSpline_beg.block(anchor_beg + 2, 0, 4, dim);
			}
			else {
				out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_beg.block(anchor_beg, 0, 4, dim);
			}
		}
		else if (segment == seg - 1) {
			out_result = T_vector.transpose() * M_matrix * G_matrixBSpline_end.block(anchor_beg, 0, 4, dim);
		}
		break;
	}
}


double cCurve::GetMaxTime() const
{
	// returns the total time needed to travel along the curve from start to end
	return mSegmentDuration * GetNumSegments();
}

void cCurve::Add(const tAnchor& anchor)
{
	mAnchors.push_back(anchor);
}

bool cCurve::ParseAnchors(const Json::Value& root)
{
	// parses an array of anchors from root

	assert(root.isArray()); 
	bool succ = true;

	int num_anchors = root.size();
	mAnchors.resize(num_anchors);

	// anchors are stored as a list of points
	// the points can be of any dimension, but the dimensions of 
	// all points should be the same
	for (int i = 0; i < num_anchors; ++i)
	{
		const auto& anchor_json = root.get(i, 0);
		tAnchor& curr_anchor = mAnchors[i];
		succ &= ParseAnchor(anchor_json, curr_anchor);
	}

	if (succ)
	{
		// compute and store the tangents at the achor points
		// these achnor tangets are currently used only for visualization
		ComputeAnchorTangents();
	}

	return succ;
}

bool cCurve::ParseAnchor(const Json::Value& root, tAnchor& out_anchor) const
{
	// parse anchors specified using a JSON format
	bool succ = true;
	if (!root[gAnchorPosKey].isNull())
	{
		const auto& pos_json = root.get(gAnchorPosKey, 0);
		int curr_dim = pos_json.size();
		out_anchor.mPos.resize(curr_dim);

		int dim = GetDim();
		succ = curr_dim == dim;
		if (!succ)
		{
			printf("Anchor dimension mismatch, expecting %i got %i\n", dim, curr_dim);
			assert(false);
		}

		if (succ)
		{
			// each anchor is defined as a list of numbers
			for (int i = 0; i < curr_dim; ++i)
			{
				out_anchor.mPos[i] = pos_json.get(i, 0).asDouble();
			}
		}
	}
	else
	{
		succ = false;
	}
	
	return succ;
}

void cCurve::PrintAnchors() const
{
	// prints out the positions of all anchors

	int num_anchors = GetNumAnchors();
	int dim = GetDim(); // dimension of each anchor
	printf("Curve Anchors:\n");
	for (int i = 0; i < num_anchors; ++i)
	{
		const tAnchor& anchor = mAnchors[i];
		printf("Anchor %i:\t", i);

		// print the position of each anchor
		for (int j = 0; j < dim; ++j)
		{
			printf("%.3f\t", anchor.mPos[j]);
		}
		printf("\n");
	}
}

void cCurve::GetAnchors(int seg, int& anchor_beg, int& anchor_end) const
{
	// compute the indices of the start and end anchors for a given curve segment
	// can be helpful when building the basis matrix
	switch (mCurveType)
	{
	case eCurveTypeCatmullRom:
		anchor_beg = seg - 1;
		anchor_end = anchor_beg + 3;
		break;
	case eCurveTypeBSpline:
		anchor_beg = seg - 2;
		anchor_end = anchor_beg + 3;
		break;
	default:
		assert(false); // unsuppoted curve type
		break;
	}
}

double cCurve::GetAnchorTime(int i) const
{
	// computes the time for a given anchor
	// i.e. roughly the time when a point will be at a particular anchor i
	int num_anchors = GetNumAnchors();
	double time = i / (num_anchors - 1.0);
	time *= GetMaxTime();
	return time;
}

void cCurve::ComputeAnchorTangents()
{
	// computes and stores the tangents at the anchor points
	for (int i = 0; i < GetNumAnchors(); ++i)
	{
		double time = GetAnchorTime(i);
		tAnchor& curr_anchor = mAnchors[i];
		EvalTangent(time, curr_anchor.mTangent);
	}
}

double cCurve::GetSegDuration(int seg) const
{
	// get the duration of each curve segment
	// for now, they are assumed to all have the same duration
	return mSegmentDuration;
}
