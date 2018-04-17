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
	
	if (value >= 1.0)
		value = 1.0;
	else if (value <= -1.0)
		value = -1.0;
	else if (isnan(value) || isinf(value))
	{
		mError(zero_vector);
		value = 1.0;
	}
	
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

//Assign x
void Node::AssignX(double x)
{
	this->coordinates(0) = x;
}

//Assign y
void Node::AssignY(double y)
{
	this->coordinates(1) = y;
}

//Assign z
void Node::AssignZ(double z)
{
	this->coordinates(2) = z;
}

//Assign id
void Node::AssignIDnumber(unsigned int i)
{
	this->IDnum = i;
}

//Assign type
void Node::AssignSubType(element_type type)
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
	double sum = 0;
	for (int i=0; i<3; i++)
		sum = sum + (this->coordinates(i) - node.coordinates(i))*(this->coordinates(i) - node.coordinates(i));
	return sqrt( sum );
}

//calc and return angle (radians)
double Node::angle(Node& n1, Node& n2)
{
	//Form vectors between (n1 and this) and (this and n2)
	Vector3D ThisN1 = *this-n1, ThisN2 = *this-n2;
	
	return ThisN1.angleRAD(ThisN2);
}

//calc and return angle (degrees)
double Node::angle_degrees(Node& n1, Node& n2)
{
	return this->angle(n1, n2) * 180.0 / M_PI;
}

//return x
double Node::getX()
{
	return this->coordinates(0);
}

//return y
double Node::getY()
{
	return this->coordinates(1);
}

//return z
double Node::getZ()
{
	return this->coordinates(2);
}

//return type
element_type Node::getType()
{
	return this->SubType;
}

//returns vector formed from subtracting nodes
Vector3D Node::operator-(const Node& node)
{
	Vector3D temp;
	for (int i=0; i<3; i++)
		temp.edit(i, this->coordinates(i) - node.coordinates(i));
	return temp;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: Node Class Definitions
 */

/*
 *								Start: LineElement Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
LineElement::LineElement()
{
	length = 0;
	IDnum = 0;
	SubType = INTERIOR;
}

//Default destructor
LineElement::~LineElement()
{
	
}

//Display info
void LineElement::DisplayInfo()
{
	std::cout << "LineElement ID: " << this->IDnum << std::endl;
	std::cout << "LineElement Type: ";
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
	std::cout << "Node 1: ( " << this->node1->getX() << " , " << this->node1->getY() << " , " << this->node1->getZ() << " )\n";
	std::cout << "Node 2: ( " << this->node2->getX() << " , " << this->node2->getY() << " , " << this->node2->getZ() << " )\n";
	std::cout << "Midpoint: ( " << this->midpoint.getX() << " , " << this->midpoint.getY() << " , " << this->midpoint.getZ() << " )\n";
	std::cout << "Length: " << this->length << std::endl << std::endl;
}

//Assign nodes
void LineElement::AssignNodes(Node& n1, Node& n2)
{
	this->node1 = &n1;
	this->node2 = &n2;
}

//Assign id
void LineElement::AssignIDnumber(unsigned int i)
{
	this->IDnum = i;
}

//Calculate length
void LineElement::calculateLength()
{
	this->length = this->node1->distance(*this->node2);
}

//Find midpoint
void LineElement::findMidpoint()
{
	this->midpoint.AssignCoordinates(0.5*(this->node1->getX()+this->node2->getX()), 0.5*(this->node1->getY()+this->node2->getY()), 0.5*(this->node1->getZ()+this->node2->getZ()));
}

//Determine type
void LineElement::determineType()
{
	if (this->node1->getType() == BOUNDARY && this->node2->getType() == BOUNDARY)
		this->midpoint.AssignSubType(BOUNDARY);
	else
		this->midpoint.AssignSubType(INTERIOR);
	this->SubType = this->midpoint.getType();
}

//eval properties
void LineElement::evaluateProperties()
{
	this->calculateLength();
	this->findMidpoint();
	this->determineType();
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: LineElement Class Definitions
 */

/*
 *								Start: SurfaceElement Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
SurfaceElement::SurfaceElement()
{
	area = 0;
	IDnum = 0;
	SubType = INTERIOR;
}

//Default destructor
SurfaceElement::~SurfaceElement()
{
	
}


//Display info
void SurfaceElement::DisplayInfo()
{
	std::cout << "SurfaceElement ID: " << this->IDnum << std::endl;
	std::cout << "SurfaceElement Type: ";
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
	std::cout << "Node 1: ( " << this->node1->getX() << " , " << this->node1->getY() << " , " << this->node1->getZ() << " )\n";
	std::cout << "Node 2: ( " << this->node2->getX() << " , " << this->node2->getY() << " , " << this->node2->getZ() << " )\n";
	std::cout << "Node 3: ( " << this->node3->getX() << " , " << this->node3->getY() << " , " << this->node3->getZ() << " )\n";
	std::cout << "Centroid: ( " << this->centroid.getX() << " , " << this->centroid.getY() << " , " << this->centroid.getZ() << " )\n";
	std::cout << "Area: " << this->area << std::endl << std::endl;
}

//Assign nodes
void SurfaceElement::AssignNodes(Node& n1, Node& n2, Node& n3)
{
	this->node1 = &n1;
	this->node2 = &n2;
	this->node3 = &n3;
}

//Assign id
void SurfaceElement::AssignIDnumber(unsigned int i)
{
	this->IDnum = i;
}

//Calculate length
void SurfaceElement::calculateArea()
{
	Vector3D n1n2 = *this->node1 - *this->node2, n2n3 = *this->node2 - *this->node3;
	Vector3D cross = n1n2.cross_product(n2n3);
	this->area = cross.norm() / 2.0;
}

//Find midpoint
void SurfaceElement::findCentroid()
{
	this->centroid.AssignCoordinates((this->node1->getX()+this->node2->getX()+this->node3->getX())/3.0, (this->node1->getY()+this->node2->getY()+this->node3->getY())/3.0, (this->node1->getZ()+this->node2->getZ()+this->node3->getZ())/3.0);
}

//Determine type
void SurfaceElement::determineType()
{
	if (this->node1->getType() == BOUNDARY && this->node2->getType() == BOUNDARY && this->node3->getType() == BOUNDARY)
		this->centroid.AssignSubType(BOUNDARY);
	else
		this->centroid.AssignSubType(INTERIOR);
	this->SubType = this->centroid.getType();
}

//eval properties
void SurfaceElement::evaluateProperties()
{
	this->calculateArea();
	this->findCentroid();
	this->determineType();
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: SurfaceElement Class Definitions
 */

/*
 *								Start: VolumeElement Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
VolumeElement::VolumeElement()
{
	volume = 0;
	IDnum = 0;
}

//Default destructor
VolumeElement::~VolumeElement()
{
	
}


//Display info
void VolumeElement::DisplayInfo()
{
	std::cout << "VolumeElement ID: " << this->IDnum << std::endl;
	std::cout << "Node 1: ( " << this->node1->getX() << " , " << this->node1->getY() << " , " << this->node1->getZ() << " )\n";
	std::cout << "Node 2: ( " << this->node2->getX() << " , " << this->node2->getY() << " , " << this->node2->getZ() << " )\n";
	std::cout << "Node 3: ( " << this->node3->getX() << " , " << this->node3->getY() << " , " << this->node3->getZ() << " )\n";
	std::cout << "Node 4: ( " << this->node4->getX() << " , " << this->node4->getY() << " , " << this->node4->getZ() << " )\n";
	std::cout << "Centroid: ( " << this->centroid.getX() << " , " << this->centroid.getY() << " , " << this->centroid.getZ() << " )\n";
	std::cout << "Volume: " << this->volume << std::endl << std::endl;
}

//Assign nodes
void VolumeElement::AssignNodes(Node& n1, Node& n2, Node& n3, Node& n4)
{
	this->node1 = &n1;
	this->node2 = &n2;
	this->node3 = &n3;
	this->node4 = &n4;
}

//Assign id
void VolumeElement::AssignIDnumber(unsigned int i)
{
	this->IDnum = i;
}

//Calculate length
void VolumeElement::calculateVolume()
{
	Vector3D n1n4 = *this->node1 - *this->node4, n2n4 = *this->node2 - *this->node4, n3n4 = *this->node3 - *this->node4;
	Vector3D cross = n2n4.cross_product(n3n4);
	this->volume = fabs(n1n4.dot_product(cross)) / 6.0;
}

//Find midpoint
void VolumeElement::findCentroid()
{
	this->centroid.AssignCoordinates((this->node1->getX()+this->node2->getX()+this->node3->getX()+this->node4->getX())/4.0, (this->node1->getY()+this->node2->getY()+this->node3->getY()+this->node4->getY())/4.0, (this->node1->getZ()+this->node2->getZ()+this->node3->getZ()+this->node4->getZ())/4.0);
}

//eval properties
void VolumeElement::evaluateProperties()
{
	this->calculateVolume();
	this->findCentroid();
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: VolumeElement Class Definitions
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
	n1.AssignIDnumber(0);
	n1.AssignCoordinates(1, 0, 0);
	n1.DisplayInfo();
	
	Node n2;
	n2.AssignSubType(INTERIOR);
	n2.AssignIDnumber(1);
	n2.AssignCoordinates(0, 1, 0);
	n2.DisplayInfo();
	
	Node n3;
	n3.AssignSubType(BOUNDARY);
	n3.AssignIDnumber(2);
	n3.AssignCoordinates(0, 0, 0);
	n3.DisplayInfo();
	
	Node n4;
	n4.AssignSubType(BOUNDARY);
	n4.AssignIDnumber(3);
	n4.AssignCoordinates(0, 0, 1);
	n4.DisplayInfo();
	
	std::cout << "Distance between n0 and n1 = " << n2.distance(n1) << std::endl;
	std::cout << "Angle (deg) between n0-n1 and n1-n2 = " << n2.angle_degrees(n1, n3) << std::endl << std::endl;;
	
	LineElement n1n2;
	n1n2.AssignNodes(n1, n2);
	n1n2.AssignIDnumber(0);
	n1n2.evaluateProperties();
	n1n2.DisplayInfo();
	
	LineElement n1n3;
	n1n3.AssignNodes(n1, n3);
	n1n3.AssignIDnumber(1);
	n1n3.evaluateProperties();
	n1n3.DisplayInfo();
	
	SurfaceElement n1n2n3;
	n1n2n3.AssignNodes(n1, n2, n3);
	n1n2n3.AssignIDnumber(0);
	n1n2n3.evaluateProperties();
	n1n2n3.DisplayInfo();
	
	VolumeElement n1n2n3n4;
	n1n2n3n4.AssignNodes(n1, n2, n3, n4);
	n1n2n3n4.AssignIDnumber(0);
	n1n2n3n4.evaluateProperties();
	n1n2n3n4.DisplayInfo();
	
	//---Exit Messages and cleanup---
	time = clock() - time;
	std::cout << "\nMESH Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	
	return success;
}
