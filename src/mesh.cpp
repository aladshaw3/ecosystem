/*!
 *  \file mesh.h
 *	\brief Mesh Objects and Associated Sub-Objects
 *	\details This kernel allows for creation of mesh objects from constitutient sub-objects
 *			such as Nodes and Elements (Line, Surface, or Volume). Mesh objects can be used
 *			to establish physical-chemcial simulations in multi-dimensional space.
 *
 *			\warning This kernel is still under active development. Use with caution!
 *
 *
 *  \author Austin Ladshaw
 *	\date 04/09/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2018, all
 *             rights reserved.
 */

#include "mesh.h"

/*
 *								Start: Vector3D Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default Constructor
Vector3D::Vector3D()
:
vector(3,1)
{
	
}

//Default Destructor
Vector3D::~Vector3D()
{
	
}

//Construction of RVector with each component
Vector3D::Vector3D(double x, double y, double z)
:
vector(3,1)
{
	vector.edit(0, 0, x);
	vector.edit(1, 0, y);
	vector.edit(2, 0, z);
}

//Copy constructor for vectors
Vector3D::Vector3D(const Vector3D &v)
:
vector(v.vector)
{
	
}

//For setting values
double& Vector3D::operator()(int i)
{
	return vector(i,0);
}

//For accessing values
double Vector3D::operator()(int i) const
{
	return vector(i,0);
}

//Compute the 2-Norm of a vector
double Vector3D::norm()
{
	return this->vector.norm();
}

//Perform the dot product between two vectors
double Vector3D::dot_product(const Vector3D &v)
{
	return this->vector.inner_product(v.vector);
}

//Returns the angle between the two vectors in radians
double Vector3D::angleRAD(Vector3D &v)
{
	double value = this->dot_product(v)/(this->norm() * v.norm());
	/*
	if (value >= 1.0)
		value = 1.0;
	else if (value <= -1.0)
		value = -1.0;
	else if (isnan(value) || isinf(value))
	{
		mError(zero_vector);
		value = 1.0;
	}
	*/
	return acos( value );
}

//Returns the angle between the two vectors in degrees
double Vector3D::angleDEG(Vector3D &v)
{
	return this->angleRAD(v) * 180.0 / M_PI;
}

//Editing a single value in the vector
void Vector3D::edit(int i, double value)
{
	this->vector.edit(i, 0, value);
}

//Editing all values in a vector
void Vector3D::set_vector(double x, double y, double z)
{
	this->vector.edit(0, 0, x);
	this->vector.edit(1, 0, y);
	this->vector.edit(2, 0, z);
}

//Vector assignment
Vector3D &Vector3D::operator=(const Vector3D &v)
{
	if (this == &v)
		return *this;
	else
	{
		this->vector=v.vector;
		
		return *this;
	}
}

//Vector Addition
Vector3D Vector3D::operator+(const Vector3D &v)
{
	Vector3D temp;
	temp.vector = this->vector + v.vector;
	return temp;
}

//Vector Subtraction
Vector3D Vector3D::operator-(const Vector3D &v)
{
	Vector3D temp;
	temp.vector = this->vector - v.vector;
	return temp;
}

//Vector dot product (short hand)
double Vector3D::operator*(const Vector3D &v)
{
	return this->dot_product(v);
}

//Vector scalar multiplication
Vector3D Vector3D::operator*(const double a)
{
	Vector3D temp;
	temp.vector = this->vector*a;
	return temp;
}

//Vector scalar division
Vector3D Vector3D::operator/(const double a)
{
	Vector3D temp;
	temp.vector = this->vector/a;
	return temp;
}

//Vector cross product
Vector3D Vector3D::cross_product(const Vector3D &v)
{
	Vector3D temp;
	temp.edit(0, this->vector(1,0)*v.vector(2,0) - this->vector(2,0)*v.vector(1,0));
	temp.edit(1, this->vector(0,0)*v.vector(2,0) - this->vector(2,0)*v.vector(0,0));
	temp.edit(2, this->vector(0,0)*v.vector(1,0) - this->vector(1,0)*v.vector(0,0));
	return temp;
}


/*
 *	-------------------------------------------------------------------------------------
 *								End: Vector3D Class Definitions
 */

/*
 *								Start: Node Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
Node::Node()
:
coordinates(0,0,0)
{
	IDnum = 0;
	SubType = INTERIOR;
	distance_tolerance = DBL_EPSILON;
}

//Default destructor
Node::~Node()
{
	
}

//Display info
void Node::DisplayInfo()
{
	std::cout << "Node ID: " << this->IDnum << std::endl;
	std::cout << "Node Type: ";
	switch (this->SubType)
	{
		case INTERIOR:
			std::cout << "INTERIOR\n";
			break;
			
		case BOUNDARY:
			std::cout << "BOUNDARY\n";
			break;
			
		default:
			std::cout << "INTERIOR\n";
			break;
	}
	std::cout << "Location: ( " << this->coordinates(0) << " , " << this->coordinates(1) << " , " << this->coordinates(2) << " )\n";
	std::cout << "Distance Tol: " << this->distance_tolerance << std::endl << std::endl;
}

//Assign coordinates
void Node::AssignCoordinates(double x, double y, double z)
{
	this->coordinates(0) = x;
	this->coordinates(1) = y;
	this->coordinates(2) = z;
}

//Assign id
void Node::AssignIDnumber(unsigned int i)
{
	this->IDnum = i;
}

//Assign type
void Node::AssignSubType(node_type type)
{
	this->SubType = type;
}

//Check for nodes in same location
bool Node::isSameLocation(Node& node)
{
	bool same = false;
	
	//Check distance
	if (this->distance(node) < this->distance_tolerance)
		same = true;
	
	return same;
}

//Check types
bool Node::isSameType(Node& node)
{
	bool same = false;
	if (this->SubType == node.SubType)
		same = true;
	return same;
}

//Check id
bool Node::isSameNode(Node& node)
{
	bool same = false;
	if (this->IDnum == node.IDnum)
		same = true;
	return same;
}

//calc and return distance
double Node::distance(Node& node)
{
	double xdiff = this->coordinates(0) - node.coordinates(0);
	double ydiff = this->coordinates(1) - node.coordinates(1);
	double zdiff = this->coordinates(2) - node.coordinates(2);
	return sqrt( xdiff*xdiff + ydiff*ydiff + zdiff*zdiff );
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: Node Class Definitions
 */

// Test function for MESH kernel
int MESH_TESTS()
{
	// --- Initializations ----
	int success = 0;
	double time = clock();
	std::cout << "\nStart Running Tests on Mesh Objects...\n\n";
	
	Node n1;
	n1.AssignSubType(BOUNDARY);
	n1.AssignIDnumber(1);
	n1.AssignCoordinates(1, 0, 0);
	n1.DisplayInfo();
	
	Node n2;
	n2.AssignSubType(INTERIOR);
	n2.AssignIDnumber(2);
	n2.AssignCoordinates(0, 0, 0);
	n2.DisplayInfo();
	
	Node n3;
	n3.AssignSubType(BOUNDARY);
	n3.AssignIDnumber(3);
	n3.AssignCoordinates(0, 1, 0);
	n3.DisplayInfo();
	
	std::cout << "Distance between n1 and n2 = " << n2.distance(n1) << std::endl;
	
	//---Exit Messages and cleanup---
	time = clock() - time;
	std::cout << "\nMESH Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	
	return success;
}
