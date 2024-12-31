/*=================================================================
 *
 * planner.cpp
 *
 *=================================================================*/
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <algorithm>

#include <tuple>
#include <string>
#include <stdexcept>
#include <regex> // For regex and split logic
#include <iostream> // cout, endl
#include <fstream> // For reading/writing files
#include <assert.h> 
#include <limits.h>
#include <unordered_map>
#include <list>
#include <queue>
#include <memory>
#include <unordered_set>

#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm
#define LINKLENGTH_CELLS 10

#ifndef MAPS_DIR
#define MAPS_DIR "../maps"
#endif
#ifndef OUTPUT_DIR
#define OUTPUT_DIR "../output"
#endif


// Some potentially helpful imports
using std::vector;
using std::array;
using std::string;
using std::runtime_error;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::cout;
using std::endl;

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                GIVEN FUNCTIONS                                                    //
//                                                                                                                   //
//*******************************************************************************************************************//

/// @brief 
/// @param filepath 
/// @return map, x_size, y_size
tuple<double*, int, int> loadMap(string filepath) {
	std::FILE *f = fopen(filepath.c_str(), "r");
	if (f) {
	}
	else {
		printf("Opening file failed! \n");
		throw runtime_error("Opening map file failed!");
	}
	int height, width;
	if (fscanf(f, "height %d\nwidth %d\n", &height, &width) != 2) {
		throw runtime_error("Invalid loadMap parsing map metadata");
	}
	
	////// Go through file and add to m_occupancy
	double* map = new double[height*width];

	double cx, cy, cz;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			char c;
			do {
				if (fscanf(f, "%c", &c) != 1) {
					throw runtime_error("Invalid parsing individual map data");
				}
			} while (isspace(c));
			if (!(c == '0')) { 
				map[y+x*width] = 1; // Note transposed from visual
			} else {
				map[y+x*width] = 0;
			}
		}
	}
	fclose(f);
	return make_tuple(map, width, height);
}

// Splits string based on deliminator
vector<string> split(const string& str, const string& delim) {   
		// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
		const std::regex ws_re(delim);
		return { std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator() };
}


double* doubleArrayFromString(string str) {
	vector<string> vals = split(str, ",");
	double* ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i) {
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double* v1, double *v2, int size) {
    for (int i = 0; i < size; ++i) {
        if (abs(v1[i]-v2[i]) > 1e-3) {
            cout << endl;
            return false;
        }
    }
    return true;
}

typedef struct {
	int X1, Y1;
	int X2, Y2;
	int Increment;
	int UsingYIndex;
	int DeltaX, DeltaY;
	int DTerm;
	int IncrE, IncrNE;
	int XIndex, YIndex;
	int Flipped;
} bresenham_param_t;

//converts float coordinates to discrete coordinates
void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size) {
	double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize)); // convert to int
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params) {
	params->UsingYIndex = 0;

	if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
		(params->UsingYIndex)++;

	if (params->UsingYIndex)
		{
			params->Y1=p1x;
			params->X1=p1y;
			params->Y2=p2x;
			params->X2=p2y;
		}
	else
		{
			params->X1=p1x;
			params->Y1=p1y;
			params->X2=p2x;
			params->Y2=p2y;
		}

	 if ((p2x - p1x) * (p2y - p1y) < 0)
		{
			params->Flipped = 1;
			params->Y1 = -params->Y1;
			params->Y2 = -params->Y2;
		}
	else
		params->Flipped = 0;

	if (params->X2 > params->X1)
		params->Increment = 1;
	else
		params->Increment = -1;

	params->DeltaX=params->X2-params->X1;
	params->DeltaY=params->Y2-params->Y1;

	params->IncrE=2*params->DeltaY*params->Increment;
	params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
	params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

	params->XIndex = params->X1;
	params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y) {
	if (params->UsingYIndex) {
        *y = params->XIndex;
        *x = params->YIndex;
        if (params->Flipped)
            *x = -*x;
    }
	else {
        *x = params->XIndex;
        *y = params->YIndex;
        if (params->Flipped)
            *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params) {
	if (params->XIndex == params->X2) {
        return 0;
    }
	params->XIndex += params->Increment;
	if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
		params->DTerm += params->IncrE;
	else {
        params->DTerm += params->IncrNE;
        params->YIndex += params->Increment;
	}
	return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
			 int x_size, int y_size) {
	bresenham_param_t params;
	int nX, nY; 
	short unsigned int nX0, nY0, nX1, nY1;

	//printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
		
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

	//printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
			return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
			 int x_size, int y_size) {
    double x0,y0,x1,y1;
    int i;
		
	 //iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
	y1 = 0;
	for(i = 0; i < numofDOFs; i++){
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
			return 0;
	}    
	return 1;
}

double* get_random_config(int numofDOFs, double* map, int x_size, int y_size) {
	double* angles = new double[numofDOFs];
	do{
		for (int i = 0; i < numofDOFs; i++) {
			//printf("generating random angle %d\n", i);
			angles[i] = ((double)rand()/(double)RAND_MAX)*2*PI;
			//printf("angle %f\n", angles[i]);
		}

	}while(!IsValidArmConfiguration(angles, numofDOFs, map, x_size, y_size));
	return angles;
}

// Forward declaration of Node
struct Node;

// Custom hash function for std::shared_ptr<Node>
struct NodePtrHash {
    std::size_t operator()(const std::shared_ptr<Node>& node) const {
        return std::hash<Node*>()(node.get());
    }
};

// Custom equality function for std::shared_ptr<Node>
struct NodePtrEqual {
    bool operator()(const std::shared_ptr<Node>& lhs, const std::shared_ptr<Node>& rhs) const {
        return lhs.get() == rhs.get();
    }
};

struct Node {

	std::shared_ptr<Node> parent;
	std::shared_ptr<Node> right;
	std::shared_ptr<Node> left;
	std::shared_ptr<Node> path_parent;
	std::vector<std::shared_ptr<Node>> path_children;
	int dof;
	double* angles;
	int id;
	int pathlength;
	double cost;
	int index;
	double distance_to_goal;

	Node(std::shared_ptr<Node> parent= nullptr, int dof=5, double* input=nullptr, std::shared_ptr<Node> right = nullptr, std::shared_ptr<Node> left = nullptr, std::shared_ptr<Node> path_parent = nullptr, int id = 0, int pathlength = 0, double cost = 0, int index = 0, double distance_to_goal = 0, std::vector<std::shared_ptr<Node>> path_children = std::vector<std::shared_ptr<Node>>()) {
		this->parent = parent;
		this->dof = dof;
		this->angles = new double[dof];
		for (int i = 0; i < dof; i++) {
			this->angles[i] = input[i];
		}
		this->right = right;
		this->left = left;
		this->path_parent = path_parent;
		this->id = id;
		this->pathlength = pathlength;
		this->cost = cost;
		this->index = index;
		this->distance_to_goal = distance_to_goal;
		this->path_children = path_children;
	}

	bool checkEqual(double* input) {
		for (int i = 0; i < dof; i++) {
			if (angles[i] != input[i]) {
				return false;
			}
		}
		return true;
	}

	double distance(double* input) {
		double dist = 0;
		double temp=0;
		for (int i = 0; i < dof; i++) {
			temp = angles[i] - input[i];
			if(fabs(temp) > PI){
				temp = 2*PI - abs(temp);
			}
			dist += pow(angles[i] - input[i], 2);
		}
		return sqrt(dist);
	}

	double* normdist(double* input){
		//calculates normalized distance between input and current node
		double* dist = new double[dof];
		double norm=0;
		for (int i = 0; i < dof; i++) {
			dist[i] = input[i] - angles[i];
			if (fabs(dist[i]) > PI){
					dist[i] = 2*PI- dist[i];
				}
			
			norm += pow(dist[i], 2);
		}
		norm = sqrt(norm);
		for (int i = 0; i < dof; i++) {
			dist[i] = dist[i]/norm;
		}

		return dist;
	}

	std::shared_ptr<Node> newConfig(double* input, double step_size, double*	map,
			 int x_size, int y_size){
		//calculates new node that is step_size away from current node in the direction of input
		double* dist = normdist(input);
		double* new_angles = new double[dof];
		double* new_angles2 = new double[dof];
		double* new_angles3 = new double[dof];
		double* new_angles4 = new double[dof];
		double* new_angles5 = new double[dof];
		for (int i = 0; i < dof; i++){
			new_angles[i] = angles[i] + step_size*dist[i];
			new_angles2[i] = angles[i] + 0.2*step_size*dist[i];
			new_angles3[i] = angles[i] + 0.4*step_size*dist[i];
			new_angles4[i] = angles[i] + 0.6*step_size*dist[i];
			new_angles5[i] = angles[i] + 0.8*step_size*dist[i];
		}
		if (IsValidArmConfiguration(new_angles, dof, map, x_size, y_size) && IsValidArmConfiguration(new_angles2, dof, map, x_size, y_size) && IsValidArmConfiguration(new_angles3, dof, map, x_size, y_size)&& IsValidArmConfiguration(new_angles4, dof, map, x_size, y_size)&& IsValidArmConfiguration(new_angles5, dof, map, x_size, y_size)){
			delete[] new_angles2;
			delete[] new_angles3;
			delete[] new_angles4;
			delete[] new_angles5;
			return std::make_shared<Node>(nullptr, dof, new_angles);
		}
		else {
			return nullptr;
		}
	}

	std::shared_ptr<Node> newConfig0(double* input, double step_size, double* map, int x_size, int y_size, int precision=4){

		double* new_angles = new double[dof];
		double* dist = new double[dof];

		dist = normdist(input);
		double factor;
		for (int k=1; k<=precision; k++){
			for (int i = 0; i < dof; i++){
				factor = (double)k/(precision+1);  // for k=2 factors: 0.33, 0.66, for k=3 factors: 0.25, 0.5, 0.75
				new_angles[i] = angles[i] + factor*dist[i];
			}
			if (!IsValidArmConfiguration(new_angles, dof, map, x_size, y_size)){
				delete[] new_angles;
				return nullptr;
			}
		}
		for (int i = 0; i < dof; i++){
			new_angles[i] = angles[i] + step_size*dist[i];
		}
		if (IsValidArmConfiguration(new_angles, dof, map, x_size, y_size)){
			delete[] new_angles;
			return std::make_shared<Node>(nullptr, dof, new_angles);
		}
		else{
			delete[] new_angles;
			return nullptr;

		}

	}

	bool checkValid(double* input, double*	map, int x_size, int y_size, int precision=5){
		double* new_angles = new double[dof];
		double* dist = new double[dof];

		if (!IsValidArmConfiguration(input, dof, map, x_size, y_size)){
			return false;
		}
		//printf("goal input valid\n");
		double factor;
		for (int k=1; k<=precision; k++){
			//printf("k %d\n", k);
			for (int i = 0; i < dof; i++){
				//printf("input %f\n", input[i]);
				//printf("angles %f\n", angles[i]);
				dist[i] = input[i] - angles[i];
				//printf("dist calculated\n");
				factor = (double)k/(precision+1);  // for k=2 factors: 0.33, 0.66, for k=3 factors: 0.25, 0.5, 0.75
				//printf("factor calculated\n");
				new_angles[i] = angles[i] + factor*dist[i];
				//printf("new angles calculated\n");
			}
			if (!IsValidArmConfiguration(new_angles, dof, map, x_size, y_size)){
				//printf("intermediate input not valid\n");
				delete[] new_angles;
				return false;
			}

		}
		// if(IsValidArmConfiguration(input, dof, map, x_size, y_size) && IsValidArmConfiguration(new_angles, dof, map, x_size, y_size) && IsValidArmConfiguration(new_angles2, dof, map, x_size, y_size)){
		// 	delete[] new_angles;
		// 	delete[] new_angles2;
		// 	return true;
		// }
		// else{
		// 	return false;
		// }
		return true;
	
	}
	

	void printNodeInfo(){
		printf("parent info: ");
		if (parent != nullptr){
			for (int i = 0; i < dof; i++){
				printf("%f ", parent->angles[i]);
			}
		}
		else{
			printf("None");
		}
		printf("\n");
		printf("current info: ");
		for (int i = 0; i < dof; i++){
			printf("%f ", angles[i]);
		}
		printf("\n");
		printf("right info: ");
		if (right != nullptr){
			for (int i = 0; i < dof; i++){
				printf("%f ", right->angles[i]);
			}
		}
		else{
			printf("None");
		}
		printf("\n");
		printf("left info: ");
		if (left != nullptr){
			for (int i = 0; i < dof; i++){
				printf("%f ", left->angles[i]);
			}
		}
		else{
			printf("None");
		}
		printf("\n");
	}

	void printNodeAngles(int level=0){
		for (int i = 0; i < dof; i++){
			printf("%f ", angles[i]);
		}
		printf("\n");
		if (left != nullptr){
			printf("left at level %d split at index %d\n", level, (id+1)%dof);
			left->printNodeAngles(level+1);
		}
		if (right != nullptr){
			printf("right at level %d split at index %d\n", level, (id+1)%dof);
			right->printNodeAngles(level+1);
		}

	}
};

struct CompareNode {
    bool operator()(const std::shared_ptr<Node>& lhs, const std::shared_ptr<Node>& rhs) const {
        return (lhs->cost + lhs->distance_to_goal) > (rhs->cost + rhs->distance_to_goal);
    }
};

// struct NodePtrHash {
//     std::size_t operator()(const std::shared_ptr<Node>& node) const {
//         return std::hash<Node*>()(node.get());
//     }
// };

// struct NodePtrEqual {
//     bool operator()(const std::shared_ptr<Node>& lhs, const std::shared_ptr<Node>& rhs) const {
//         return lhs.get() == rhs.get();
//     }
// };

struct Tree {

	std::shared_ptr<Node> root;
	int dof;
	double* angles;

	/**
	 * @brief Constructs a Tree object with a given root node, degrees of freedom (dof), and input angles.
	 * 
	 * @param root A shared pointer to the root node of the tree.
	 * @param dof The number of degrees of freedom.
	 * @param input A pointer to an array of input angles.
	 */
	Tree(std::shared_ptr<Node> root, int dof, double* input) {
		this->root = root;
		this->dof = dof;
		this->angles = new double[dof];
		for (int i = 0; i < dof; i++) {
			this->angles[i] = input[i];
		}
	}

	void printTree(){
		root->printNodeAngles();
	}

	void printNumNodes(){
		int count = 0;
		std::queue<std::shared_ptr<Node>> q;
		q.push(root);
		while (!q.empty()) {
			std::shared_ptr<Node> curr_node = q.front();
			q.pop();
			count++;
			if (curr_node->left != nullptr) {
				q.push(curr_node->left);
			}
			if (curr_node->right != nullptr) {
				q.push(curr_node->right);
			}
		}
		printf("-----------------Number of nodes in tree: %d\n", count);
	}

	/**
	 * @brief Inserts a new node into the binary tree.
	 *
	 * This function inserts a new node into the binary tree based on the angles
	 * of the nodes. The tree is traversed starting from the root, and the new node
	 * is inserted either to the left or right of the current node depending on the
	 * comparison of their angles. The function ensures that the new node is inserted
	 * at the correct position and updates its parent and id accordingly.
	 *
	 * @param newNode A shared pointer to the new node to be inserted.
	 */
	void insert(std::shared_ptr<Node> newNode){

		std::shared_ptr<Node> curr_node = root;
		int index = 0;
		while(curr_node!=nullptr){
			if (newNode->angles[index] < curr_node->angles[index]){
				if (curr_node->left == nullptr){
					curr_node->left = newNode;
					newNode->id = (index+1)%dof;
					newNode->parent = curr_node;
					//printf("new node at index %d %f\n", newNode->id, newNode->angles[index]);
					break;
				}else{
					curr_node = curr_node->left;
				}
			}else{
				if (curr_node->right == nullptr){
					curr_node->right = newNode;
					newNode->id = (index+1)%dof;
					newNode->parent = curr_node;
					//printf("new node at index %d %f\n", newNode->id, newNode->angles[index]);
					break;
				}else{
					curr_node = curr_node->right;
				}
			}
			index = (index+1)%dof;

		}

	}

	/**
	 * @brief Searches for a node with the specified input values in the tree.
	 * 
	 * This function traverses the binary tree starting from the root node and 
	 * checks if there exists a node with the same values as the input array.
	 * 
	 * @param input A pointer to an array of double values representing the input to search for.
	 * @return true if a node with the specified input values is found, false otherwise.
	 */
	bool find(double* input){
		std::shared_ptr<Node> curr_node = root;
		while(curr_node != nullptr){
			if (curr_node->checkEqual(input)){
				return true;
			}
			if (input[curr_node->id] < curr_node->angles[curr_node->id]){
				curr_node = curr_node->left;
			}else{
				curr_node = curr_node->right;
			}
		}
		return false;
	}

	/**
	 * @brief Finds a node in the tree that matches the given input.
	 *
	 * This function traverses the tree starting from the root node and searches for a node
	 * whose angles match the given input. If a matching node is found, it is returned.
	 * Otherwise, the function returns nullptr.
	 *
	 * @param input A pointer to an array of doubles representing the input angles.
	 * @return A shared pointer to the found Node if a match is found, otherwise nullptr.
	 */
	std::shared_ptr<Node> findNode(double* input){
		std::shared_ptr<Node> curr_node = root;
		while(curr_node != nullptr){
			if (curr_node->checkEqual(input)){
				return curr_node;
			}
			if (input[curr_node->id] < curr_node->angles[curr_node->id]){
				curr_node = curr_node->left;
			}else{
				curr_node = curr_node->right;
			}
		}
		return nullptr;
	}

	/**
	 * @brief Computes the depth of the given node in the tree.
	 *
	 * This function traverses up the tree from the given node to the root,
	 * counting the number of edges (or levels) between the given node and the root.
	 *
	 * @param curr_node A shared pointer to the current node whose depth is to be calculated.
	 * @return The depth of the given node in the tree.
	 */
	int getDepth(std::shared_ptr<Node> curr_node){
		int depth = 0;
		while(curr_node->parent != nullptr){
			curr_node = curr_node->parent;
			depth++;
		}
		return depth;
	}

	/**
	 * @brief Finds the nearest neighbor node to the given input in a k-d tree.
	 * 
	 * This function recursively searches through a k-d tree to find the node that is closest
	 * to the given input point. It updates the minimum distance found so far and returns
	 * the nearest neighbor node.
	 * 
	 * @param input A pointer to an array of doubles representing the input point.
	 * @param curr_node A shared pointer to the current node in the k-d tree.
	 * @param mind A reference to a double representing the minimum distance found so far.
	 * @return A shared pointer to the nearest neighbor node.
	 */
	std::shared_ptr<Node> nearestNeighbor(double* input, std::shared_ptr<Node> curr_node, double& mind){

		std::shared_ptr<Node> NN;
		if (curr_node != nullptr){

			double curr_dist = curr_node->distance(input);
			
			if (curr_dist < mind) {
				mind = curr_dist;
				NN = curr_node;
			}else{
				NN = nullptr;
			}
			int index = curr_node->id;
			std::shared_ptr<Node> left = curr_node->left;
			std::shared_ptr<Node> right = curr_node->right;

			if (input[index] < curr_node->angles[index]){
				//printf("going left\n");
				if(left != nullptr){
					double temp = mind;
					//printf("calling again\n");	
					std::shared_ptr<Node> leftNN = nearestNeighbor(input, left, temp);
					if (leftNN != nullptr){
						mind = temp;
						NN = leftNN;
					}
				}
				if (right != nullptr){
					if (right->distance(input) < mind){
						//printf("calling again\n");
						std::shared_ptr<Node> rightNN = nearestNeighbor(input, right, mind);
						NN = rightNN;
					}
				}

				
			} else {
				//printf("going right\n");
				if (right!=nullptr){
					double temp = mind;
					//printf("calling again\n");
					std::shared_ptr<Node> rightNN = nearestNeighbor(input, right, temp);
					if (rightNN != nullptr){
						mind = temp;
						NN = rightNN;
					}
				}
				if (left != nullptr){
					if (left->distance(input) < mind){
						//printf("calling again\n");
						std::shared_ptr<Node> leftNN = nearestNeighbor(input, left, mind);
						NN = leftNN;
					}
				}
			}

		}
		return NN;
	}

	/**
	 * @brief Finds the nearest neighbor node to the given input in a tree structure.
	 *
	 * This function recursively searches for the nearest neighbor node to the given input
	 * starting from the current node. It updates the minimum distance found so far and 
	 * returns the nearest neighbor node.
	 *
	 * @param input A pointer to an array of doubles representing the input coordinates.
	 * @param curr_node A shared pointer to the current node in the tree.
	 * @param mind A reference to a double representing the minimum distance found so far.
	 * @return A shared pointer to the nearest neighbor node. If no nearest neighbor is found, returns nullptr.
	 */
	std::shared_ptr<Node> nearestNeighbor2(double* input, std::shared_ptr<Node> curr_node, double& mind){

		std::shared_ptr<Node> NN;
		//printf("NearestNeighbor\n");
		double curr_dist;
		if (curr_node != nullptr){
			if (!curr_node->checkEqual(input)){
				curr_dist = curr_node->distance(input);
				if (curr_dist == 0) {
					//printf("found equal\n");
					NN = nullptr;}
				else{ if (curr_dist < mind) {
					//printf("found closer\n");
					//printf("curr dist %f\n", curr_dist);
					mind = curr_dist;
					NN = curr_node;
				}else{
					//printf("found worse\n");
					NN = nullptr;
				}}
			}
			else{
				NN = nullptr;
				//printf("found equal\n");
			}
			
			int index = curr_node->id;
			std::shared_ptr<Node> left = curr_node->left;
			std::shared_ptr<Node> right = curr_node->right;

			if (input[index] < curr_node->angles[index]){
				//printf("going left\n");
				if(left != nullptr){
					double temp = mind;
					//printf("calling again\n");	
					std::shared_ptr<Node> leftNN = nearestNeighbor2(input, left, temp);
					if (leftNN != nullptr){
						mind = temp;
						NN = leftNN;
					}
				}
				if (right != nullptr){
					if (right->distance(input) < mind){
						//printf("calling again\n");
						std::shared_ptr<Node> rightNN = nearestNeighbor2(input, right, mind);
						NN = rightNN;
					}
				}

				
			} else {
				//printf("going right\n");
				if (right!=nullptr){
					double temp = mind;
					//printf("calling again\n");
					std::shared_ptr<Node> rightNN = nearestNeighbor2(input, right, temp);
					if (rightNN != nullptr){
						mind = temp;
						NN = rightNN;
					}
				}
				if (left != nullptr){
					if (left->distance(input) < mind){
						//printf("calling again\n");
						std::shared_ptr<Node> leftNN = nearestNeighbor2(input, left, mind);
						NN = leftNN;
					}
				}
			}
		}
		return NN;
	}

	/**
	 * @brief Finds all nodes within a specified range of the input point and adds them to the neighbors vector.
	 * 
	 * This function recursively traverses the tree to find all nodes that are within the given range of the input point.
	 * If a node is within the range and is not equal to the input point, it is added to the neighbors vector.
	 * 
	 * @param input A pointer to the input point coordinates.
	 * @param curr_node The current node being checked.
	 * @param range The distance range within which neighbors are to be found.
	 * @param neighbors A vector to store the pointers to the neighboring nodes found within the range.
	 */
	void nearestRangeNeighbor(double* input, std::shared_ptr<Node> curr_node, double range, std::vector<std::shared_ptr<Node>>& neighbors){
		
		if (curr_node == nullptr){
			return;
		}
		double curr_dist = curr_node->distance(input);
		if (curr_dist < range){
			if (!curr_node->checkEqual(input)){
				neighbors.push_back(curr_node);
			}
			//neighbors.push_back(curr_node);
		}
		int index = curr_node->id;

		if (input[index] < curr_node->angles[index]){
			nearestRangeNeighbor(input, curr_node->left, range, neighbors);
			if (curr_node->angles[index] - input[index] < range){
				nearestRangeNeighbor(input, curr_node->right, range, neighbors);
			}
		} else {
			nearestRangeNeighbor(input, curr_node->right, range, neighbors);
			if (input[index] - curr_node->angles[index] < range){
				nearestRangeNeighbor(input, curr_node->left, range, neighbors);
			}
		}

	}

    void extend(double* goal, double* qrand, double step_size, double* map, int x_size, int y_size){
		
		double mind = INT_MAX;
		std::shared_ptr<Node> NN = nearestNeighbor(qrand, root, mind);
		if ((NN->distance(qrand) < step_size)&& (NN->checkValid(qrand, map, x_size, y_size))){
			std::shared_ptr<Node> NewN=std::make_shared<Node>(nullptr, dof, qrand);
			NewN->path_parent = NN;
			NewN->pathlength = NN->pathlength + 1;
			NewN->cost = NN->cost + NewN->distance(NN->angles);
			insert(NewN);
			return;
		}
		//printf("found nearest neighbor\n");
		std::shared_ptr<Node> q_near = NN->newConfig(qrand, step_size, map, x_size, y_size);
		//printf("found new config\n");
		if (q_near != nullptr){
			q_near->path_parent = NN;
			q_near->pathlength = NN->pathlength + 1;
			q_near->cost = NN->cost + q_near->distance(NN->angles);
			insert(q_near);
			if ((q_near->distance(goal) < step_size)&& (q_near->checkValid(goal, map, x_size, y_size))){
				std::shared_ptr<Node> NewN=std::make_shared<Node>(nullptr, dof, goal);
				NewN->path_parent = NN;
				NewN->pathlength = NN->pathlength + 1;
				NewN->cost = NN->cost + NewN->distance(NN->angles);
				insert(NewN);
				
			}
			else{
				//printf("advanced\n");
			}
		}
		else{
			//printf("trapped\n");
		}

	}

	/**
	 * @brief Extends the path towards a random configuration or goal configuration.
	 *
	 * This function attempts to extend the path from the nearest neighbor node towards
	 * a random configuration (qrand) or the goal configuration. If the random configuration
	 * is within the step size and is valid, it is added to the path. Otherwise, a new 
	 * configuration is generated towards the random configuration within the step size.
	 * If this new configuration is valid and close enough to the goal, the goal is added
	 * to the path.
	 *
	 * @param goal The goal configuration.
	 * @param qrand The random configuration.
	 * @param step_size The maximum step size for extending the path.
	 * @param map The map of the environment.
	 * @param x_size The size of the map in the x dimension.
	 * @param y_size The size of the map in the y dimension.
	 * @return A pointer to the new configuration if the path is extended, otherwise nullptr.
	 */
	double* extend2(double* goal, double* qrand, double step_size, double* map, int x_size, int y_size){
		
		double mind = INT_MAX;
		double *qnew;
		std::shared_ptr<Node> NN = nearestNeighbor(qrand, root, mind);
		if ((NN->distance(qrand) < step_size)&& (NN->checkValid(qrand, map, x_size, y_size))){
			std::shared_ptr<Node> NewN=std::make_shared<Node>(nullptr, dof, qrand);
			NewN->path_parent = NN;
			NewN->pathlength = NN->pathlength + 1;
			NewN->cost = NN->cost + NewN->distance(NN->angles);
			insert(NewN);
			return qrand;
		}
		//printf("found nearest neighbor\n");
		std::shared_ptr<Node> q_near = NN->newConfig(qrand, step_size, map, x_size, y_size);
		//printf("found new config\n");
		if (q_near != nullptr){
			q_near->path_parent = NN;
			q_near->pathlength = NN->pathlength + 1;
			q_near->cost = NN->cost + q_near->distance(NN->angles);
			insert(q_near);
			if ((q_near->distance(goal) < step_size)&& (q_near->checkValid(goal, map, x_size, y_size))){
				std::shared_ptr<Node> NewN=std::make_shared<Node>(nullptr, dof, goal);
				NewN->path_parent = NN;
				NewN->pathlength = NN->pathlength + 1;
				NewN->cost = NN->cost + NewN->distance(NN->angles);
				insert(NewN);
				qnew = new double[dof];
				qnew=goal;
				//printf("Reached goal\n");
			}
			else{
				//printf("advanced\n");
				qnew = new double[dof];
				qnew = q_near->angles;
			}
		}
		else{
			//printf("trapped\n");
			qnew=nullptr;
		}
		return qnew;

	}

	/**
	 * @brief Attempts to extend the path towards a goal configuration.
	 *
	 * This function tries to extend the path from the nearest node towards a random configuration (qrand) 
	 * by a given step size. If the extension reaches the goal configuration, it inserts the goal node 
	 * into the tree. If the extension is valid but does not reach the goal, it inserts the new configuration 
	 * into the tree and returns the new node. If the extension is not valid, it returns nullptr.
	 *
	 * @param goal Pointer to the goal configuration.
	 * @param qrand Pointer to the random configuration.
	 * @param step_size The step size for the extension.
	 * @param map Pointer to the map data.
	 * @param x_size The size of the map in the x dimension.
	 * @param y_size The size of the map in the y dimension.
	 * @param Nnode Optional parameter specifying a node to start the extension from. If nullptr, the nearest 
	 *              neighbor to qrand is used.
	 * @return A shared pointer to the nearest node reached by the extension, or nullptr if the extension is not valid.
	 */
	std::shared_ptr<Node> connect_extend(double* goal, double* qrand, double step_size, double* map, int x_size, int y_size, std::shared_ptr<Node> Nnode=nullptr){
		
		double mind = INT_MAX;
		std::shared_ptr<Node> nearest;
		std::shared_ptr<Node> NN;
		if (Nnode == nullptr){
			NN = nearestNeighbor(qrand, root, mind);
		}
		else{
			NN = Nnode;
		}

		if ((NN->distance(qrand) < step_size)&& (NN->checkValid(qrand, map, x_size, y_size))){
			std::shared_ptr<Node> goalnode=std::make_shared<Node>(nullptr, dof, qrand);
			goalnode->path_parent = NN;
			goalnode->pathlength = NN->pathlength + 1;
			insert(goalnode);
			nearest=goalnode; //reached
			return nearest;
		}
		//printf("found nearest neighbor\n");
		std::shared_ptr<Node> q_near = NN->newConfig(qrand, step_size, map, x_size, y_size);
		//printf("found new config\n");
		if (q_near != nullptr){
			q_near->path_parent = NN;
			q_near->pathlength = NN->pathlength + 1;
			insert(q_near);
			if ((q_near->distance(goal) < step_size)&& (q_near->checkValid(goal, map, x_size, y_size) && !find(goal))){
				std::shared_ptr<Node> goalnode=std::make_shared<Node>(nullptr, dof, goal);
				goalnode->path_parent = q_near;
				goalnode->pathlength = q_near->pathlength + 1;
				insert(goalnode);
				nearest=q_near;//goalnode; //maybe?
				printf("-----Reached goal\n");
			}
			else{
				nearest= q_near;
				//printf("advanced\n");
			}
		}
		else{
			nearest = nullptr;
			//printf("trapped\n");
		}
		
		return nearest;
	}

	/**
	 * @brief Attempts to connect to a goal configuration from a random configuration using a step size.
	 *
	 * This function tries to connect to the goal configuration from a given random configuration (qrand)
	 * by iteratively extending towards the goal using a specified step size. It uses a map to check for
	 * valid configurations and updates the nearest node in each iteration. The process continues until
	 * the goal is reached or no further progress can be made.
	 *
	 * @param goal Pointer to the goal configuration.
	 * @param qrand Pointer to the random configuration.
	 * @param step_size The step size used for extending towards the goal.
	 * @param map Pointer to the map used for checking valid configurations.
	 * @param x_size The size of the map in the x dimension.
	 * @param y_size The size of the map in the y dimension.
	 * @return True if the goal configuration is reached, false otherwise.
	 */
	bool connect(double* goal, double* qrand, double step_size, double* map, int x_size, int y_size){
		bool reached=false;
		std::shared_ptr<Node> nearest = nullptr;
		int count = 0;
		do{
			nearest= connect_extend(goal, qrand, step_size, map, x_size, y_size, nearest);
			//printf("found nearest\n");
			//printf("nearest angles %f %f %f %f %f\n", nearest->angles[0], nearest->angles[1], nearest->angles[2], nearest->angles[3], nearest->angles[4]);
			count++;
			if (count%1000==0){
				printf("count %d\n", count);
				printf("nearest angles %f %f %f %f %f\n", nearest->angles[0], nearest->angles[1], nearest->angles[2], nearest->angles[3], nearest->angles[4]);
				printf("qrand angles %f %f %f %f %f\n", qrand[0], qrand[1], qrand[2], qrand[3], qrand[4]);

			}

			if (nearest!=nullptr){
				if (nearest->checkEqual(qrand)){
					break;
				}
			}
		} while (nearest!=nullptr);

		if (nearest != nullptr){
			
			if (nearest->checkEqual(qrand)){
				printf("reached\n");
				reached=true;
			}
		}
		else {
			//printf("trapped\n");
		}

		//maybe check for goal? idk
		return reached;
		
		}

	/**
	 * @brief Builds a Rapidly-exploring Random Tree (RRT) to find a path to the goal.
	 * 
	 * This function attempts to build an RRT by sampling random configurations and extending
	 * the tree towards the goal. It prints progress information and stops when the goal is found
	 * or the maximum number of samples is reached.
	 * 
	 * @param goal A pointer to the goal configuration.
	 * @param num_samples The maximum number of samples to generate.
	 * @param step_size The step size for extending the tree.
	 * @param map A pointer to the map data.
	 * @param x_size The width of the map.
	 * @param y_size The height of the map.
	 */

	void build_RRT(double* goal, int num_samples, double step_size, double* map, int x_size, int y_size) {
		int k = 0;
		double mind;
		printf("building RRT\n");
		while (k<num_samples){
			double* qrand = get_random_config(dof, map, x_size, y_size);
			mind = INT_MAX;
			if (k%500000==0){
				printf("---------sample %d\n", k);
				std::shared_ptr<Node> closest_togoal=nearestNeighbor(goal, root, mind);
				if (closest_togoal!=nullptr){
					printf("closest distance to goal %f\n", closest_togoal->distance(goal));
				}
				else{
					printf("closest to goal is null\n");
				}
			}
			mind = INT_MAX;
			extend(goal, qrand, step_size, map, x_size, y_size);
			k++;
			if(find(goal)){
				printf("goal found after %d samples\n", k);
				break;
			}

		}
	}

	/**
	 * @brief Implements the RRT* (Rapidly-exploring Random Tree Star) algorithm for path planning.
	 * 
	 * This function attempts to find a path from the root node to the goal configuration using the RRT* algorithm.
	 * 
	 * @param goal A pointer to an array representing the goal configuration.
	 * @param num_samples The number of samples to be generated.
	 * @param step_size The step size to be used for extending the tree.
	 * @param map A pointer to the map data.
	 * @param x_size The width of the map.
	 * @param y_size The height of the map.
	 * 
	 * The function generates random configurations and attempts to connect them to the nearest node in the tree.
	 * It uses a goal bias to occasionally sample the goal configuration directly.
	 * The tree is extended towards the random configuration or the goal configuration if it is within the step size.
	 * The function also re-wires the tree to ensure that the path cost is minimized.
	 * The process continues until the goal is found or the maximum number of samples is reached.
	 */
	void RRTstar(double* goal, int num_samples, double step_size, double* map, int x_size, int y_size) {
		int k = 0;
		printf("building RRT******\n");
		double mind;
		std::shared_ptr<Node> NN;
		std::shared_ptr<Node> q_new;
		double goalbias=0.13;
		double closest_dist=INT_MAX;
		while(k<num_samples){
			double* qrand = get_random_config(dof, map, x_size, y_size);
			if((double)rand() / RAND_MAX < goalbias){
				//printf("goal bias\n");
				qrand = new double[dof];
				for(int i = 0; i < dof; i++){
					qrand[i] = goal[i];
				}
			}
			
			mind = INT_MAX;
			if (k%7500==0){
				printf("---------sample %d\n", k);
				std::shared_ptr<Node> closest_togoal=nearestNeighbor(goal, root, mind);
				if (closest_togoal!=nullptr){
					closest_dist = closest_togoal->distance(goal);
					if ((closest_dist < 8*step_size)&&(goalbias < 0.2)){
						goalbias = 0.2;
					}
					if ((closest_dist < 6*step_size)&&(goalbias < 0.3)){
						goalbias = 0.3;
					}
					printf("closest distance to goal %f\n", closest_dist);
				}
				else{
					printf("closest to goal is null\n");
				}
			}
			mind = INT_MAX;
			NN = nearestNeighbor(qrand, root, mind);
			if (NN->distance(qrand) < step_size){
				q_new = std::make_shared<Node>(nullptr, dof, qrand);
				if (!NN->checkValid(qrand, map, x_size, y_size, 5)){
					q_new= nullptr;
				}
			}
			else{
				q_new = NN->newConfig(qrand, step_size, map, x_size, y_size);
			}
			if(q_new==nullptr){
				k++;
				continue;
			}
			q_new->path_parent = NN;
			q_new->pathlength = NN->pathlength + 1;
			q_new->cost = NN->cost + q_new->distance(NN->angles);
			insert(q_new);
			if ((q_new->distance(goal) < step_size) && (q_new->checkValid(goal, map, x_size, y_size, 6))){
				std::shared_ptr<Node> NewN=std::make_shared<Node>(nullptr, dof, goal);
				NewN->path_parent = NN;
				NewN->pathlength = NN->pathlength + 1;
				NewN->cost = NN->cost + NewN->distance(NN->angles);
				insert(NewN);
				//insert(std::make_shared<Node>(nullptr, dof, goal));
				//printf("Reached goal\n");
			}
			//printf("qnew found\n");
			std::vector<std::shared_ptr<Node>> neighbors;
			nearestRangeNeighbor(q_new->angles, root, 1.96*step_size, neighbors); // 3-> 2
			//printf("neighbors found\n");
			std::shared_ptr<Node> min_node = NN;
			for (int i = 0; i < neighbors.size(); i++){
				//printf("neighbor %d\n", i);
				if (neighbors[i]->checkValid(q_new->angles, map, x_size, y_size, 5)){
					//printf("valid\n");
					if (neighbors[i]->cost + neighbors[i]->distance(q_new->angles) < q_new->cost){
						min_node = neighbors[i];
						q_new->path_parent = min_node;
						q_new->pathlength = min_node->pathlength + 1;
						q_new->cost = min_node->cost + q_new->distance(min_node->angles); //Maybe add intermediate nodes
						min_node->path_children.push_back(q_new);
					}
				}
			}
			
			for (int i = 0; i < neighbors.size(); i++){

				if (neighbors[i]->checkEqual(min_node->angles)){
					continue;
				}

				if (neighbors[i]->checkValid(q_new->angles, map, x_size, y_size, 5)){
					if (q_new->cost + q_new->distance(neighbors[i]->angles) < neighbors[i]->cost){
						std::shared_ptr<Node> old_parent = neighbors[i]->path_parent;
						if(old_parent != nullptr){
							old_parent->path_children.erase(std::remove(old_parent->path_children.begin(), old_parent->path_children.end(), neighbors[i]), old_parent->path_children.end());
						}
						neighbors[i]->path_parent = q_new;
						neighbors[i]->pathlength = q_new->pathlength + 1;
						neighbors[i]->cost = q_new->cost + neighbors[i]->distance(q_new->angles);
						q_new->path_children.push_back(neighbors[i]);
						std::shared_ptr<Node> temp = neighbors[i];
						std::vector<std::shared_ptr<Node>> temp_children = temp->path_children;
						for (int j = 0; j < temp_children.size(); j++){
							temp_children[j]->cost = temp->cost + temp_children[j]->distance(temp->angles);
						}
						
					}
				}
			}
			if(find(goal)){
				printf("goal found after %d samples\n", k+1);
				break;
			}
			k++;

		}
	}


	/**
	 * @brief Generates a plan from the goal node to the start node.
	 *
	 * This function constructs a plan by tracing back from the goal node to the start node,
	 * storing the angles at each node in the plan array. The plan is allocated dynamically
	 * and must be freed by the caller.
	 *
	 * @param goal A shared pointer to the goal node from which the plan is generated.
	 * @param plan A pointer to a double pointer array where the plan will be stored.
	 *             The plan array will be allocated within this function.
	 * @param depth The depth of the plan, indicating the number of steps from the start node to the goal node.
	 *
	 * @throws std::runtime_error If memory allocation for any step in the plan fails.
	 */
	void get_plan(std::shared_ptr<Node> goal, double*** plan, int depth){
		*plan = (double**) malloc((depth+1)*sizeof(double*));
		std::shared_ptr<Node> curr_node = goal;
		int i;
		for (int k = 0; k < depth+1; k++){
			i=depth-k;
			//printf("curr node angles %f %f %f %f %f\n", curr_node->angles[0], curr_node->angles[1], curr_node->angles[2], curr_node->angles[3], curr_node->angles[4]);

			if (curr_node == nullptr){
				printf("curr node is null\n");
				break;
			}

			(*plan)[i] = (double*) malloc(dof*sizeof(double));
			if ((*plan)[i] == nullptr) {
				throw std::runtime_error("Failed to allocate memory for plan step");
			}

			for (int j = 0; j < dof; j++){
				(*plan)[i][j] = curr_node->angles[j];
			}
			curr_node = curr_node->path_parent;	
		}
		
	}

	/**
	 * @brief Generates a plan by traversing from the goal node to the root node.
	 *
	 * This function fills the provided plan with the angles from each node in the path
	 * from the goal node to the root node. The plan is a 2D array where each row represents
	 * the angles at a specific depth in the path.
	 *
	 * @param goal A shared pointer to the goal node from which the plan is generated.
	 * @param plan A pointer to a 2D array that will be filled with the angles from each node.
	 * @param depth The depth of the path from the root node to the goal node.
	 *
	 * @throws std::runtime_error If memory allocation for the plan step fails.
	 */
	void get_plan_down(std::shared_ptr<Node> goal, double*** plan, int depth){
		std::shared_ptr<Node> curr_node = goal;
		int i;
		for (int k = 0; k < depth+1; k++){
			i=depth-k;
			//printf("curr node angles %f %f %f %f %f\n", curr_node->angles[0], curr_node->angles[1], curr_node->angles[2], curr_node->angles[3], curr_node->angles[4]);

			if (curr_node == nullptr){
				printf("curr node is null\n");
				break;
			}

			(*plan)[i] = (double*) malloc(dof*sizeof(double));
			if ((*plan)[i] == nullptr) {
				throw std::runtime_error("Failed to allocate memory for plan step");
			}

			for (int j = 0; j < dof; j++){
				(*plan)[i][j] = curr_node->angles[j];
			}
			curr_node = curr_node->path_parent; //parent;	
		}
		
	}

	/**
	 * @brief Generates a plan from the goal node up to a specified depth.
	 *
	 * This function traces back from the goal node to its ancestors up to a specified depth
	 * and stores the angles of each node in the provided plan array.
	 *
	 * @param goal A shared pointer to the goal node from which the plan is generated.
	 * @param plan A pointer to a 3D array where the plan will be stored. The array should be pre-allocated.
	 * @param depth1 The starting depth index for the plan.
	 * @param depth2 The ending depth index for the plan.
	 *
	 * @throws std::runtime_error If memory allocation for the plan step fails.
	 */
	void get_plan_up(std::shared_ptr<Node> goal, double*** plan, int depth1, int depth2){

		if (goal->parent == nullptr){
			printf("goal parent is null\n");
			return;
		}

		std::shared_ptr<Node> curr_node =  goal->path_parent; //parent
		//int i;
		for (int i = depth1+1; i < depth2; i++){
			//printf("curr node angles %f %f %f %f %f\n", curr_node->angles[0], curr_node->angles[1], curr_node->angles[2], curr_node->angles[3], curr_node->angles[4]);

			if (curr_node == nullptr){
				printf("curr node is null\n");
				break;
			}

			(*plan)[i] = (double*) malloc(dof*sizeof(double));
			if ((*plan)[i] == nullptr) {
				throw std::runtime_error("Failed to allocate memory for plan step");
			}

			for (int j = 0; j < dof; j++){
				(*plan)[i][j] = curr_node->angles[j];
			}
			curr_node = curr_node->path_parent; //parent;	
		}
		
	}

};


/**
 * @brief Implements the RRT-Connect algorithm to find a path between two trees.
 *
 * This function attempts to connect two trees, A and B, to find a path from the start to the goal configuration.
 * It iteratively samples random configurations and extends the trees towards each other until a path is found or
 * the maximum number of samples is reached.
 *
 * @param A Pointer to the first tree.
 * @param B Pointer to the second tree.
 * @param goal Pointer to the goal configuration.
 * @param num_samples Maximum number of samples to attempt.
 * @param step_size Step size for extending the trees.
 * @param map Pointer to the map data.
 * @param x_size Size of the map in the x dimension.
 * @param y_size Size of the map in the y dimension.
 * @return Pointer to the configuration where the trees connect if a path is found, otherwise nullptr.
 */
double* RRTconnect(Tree* A, Tree* B,double* goal, int num_samples, double step_size, double* map, int x_size, int y_size){
		
		Tree* currTree = A;
		Tree* otherTree=B;
		Tree* temp;
		bool reached=false;
		int k = 0;
		while (k<num_samples){
			//printf("sample %d\n", k);
			double* qrand = get_random_config(currTree->dof, map, x_size, y_size);
			//printf("got random config\n");
			double* qnear = currTree->extend2(goal, qrand, step_size, map, x_size, y_size);
			//printf("extended\n");
			if (qnear != nullptr){
				reached = otherTree->connect(goal, qnear, step_size, map, x_size, y_size);
				//printf("connected\n");
				if(reached){
					printf("reached and common path found after %d samples\n", k+1);
					return qnear;
					break;
				}
				temp = currTree;
				currTree=otherTree;
				otherTree=temp;
			
			}
			k++;
		}
		printf("No path found done\n");
		return nullptr;
}

void RRTconnectPath(Tree* A, Tree* B, double* link, double step_size, double* map, int x_size, int y_size, double*** plan, int* planlength){

	std::shared_ptr<Node> goalA = A->findNode(link);
	std::shared_ptr<Node> goalB = B->findNode(link);
	int depthA = goalA->pathlength; //A->getDepth(goalA);
	int depthB = goalB->pathlength; //B->getDepth(goalB);
	int total_depth = depthA + depthB+2;
	*planlength = total_depth-1;
	*plan = (double**) malloc(total_depth*sizeof(double*));

	A->get_plan_down(goalA, plan, depthA);
	B->get_plan_up(goalB, plan, depthA, total_depth-1);

		
}

/**
 * @brief A structure representing a graph consisting of nodes and edges.
 */
struct Graph {
	std::vector<std::shared_ptr<Node>> nodes; ///< A vector of nodes in the graph.
	std::unordered_map<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, NodePtrHash, NodePtrEqual> edges; ///< An adjacency list representing edges between nodes.

	/**
	 * @brief Default constructor for the Graph structure.
	 */
	Graph() {}

	/**
	 * @brief Adds a node to the graph.
	 * @param node A shared pointer to the node to be added.
	 */
	void addNode(std::shared_ptr<Node> node) {
		nodes.push_back(node);
	}

	/**
	 * @brief Adds an undirected edge between two nodes in the graph.
	 * @param from A shared pointer to the starting node.
	 * @param to A shared pointer to the ending node.
	 */
	void addEdge(std::shared_ptr<Node> from, std::shared_ptr<Node> to) {
		//check if from and to are already in the graph
		
		if (edges.find(from) == edges.end()) {
			edges[from] = std::vector<std::shared_ptr<Node>>();
		}
		if (std::find(edges[from].begin(), edges[from].end(), to) == edges[from].end()) {
			edges[from].push_back(to);
		}
		//edges[from].push_back(to);

		if (edges.find(to) == edges.end()) {
			edges[to] = std::vector<std::shared_ptr<Node>>();
		}
		if (std::find(edges[to].begin(), edges[to].end(), from) == edges[to].end()) {
			edges[to].push_back(from);
		}
		//edges[to].push_back(from);
	}

	/**
	 * @brief Retrieves the neighbors of a given node.
	 * @param node A shared pointer to the node whose neighbors are to be retrieved.
	 * @return A vector of shared pointers to the neighboring nodes.
	 */
	std::vector<std::shared_ptr<Node>> getNeighbors(std::shared_ptr<Node> node) {
		if (edges.find(node) != edges.end()) {
			return edges[node];
		}
		printf("no neighbors found\n");
		return std::vector<std::shared_ptr<Node>>();
	}

	/**
	 * @brief Checks if a node with the given input exists in the graph.
	 * @param input A pointer to the input data to be checked.
	 * @return True if the node exists, false otherwise.
	 */
	bool inGraph(double* input){ 
		//return edges.find(node) != edges.end();
		for (int i = 0; i < nodes.size(); i++){
			if (nodes[i]->checkEqual(input)){
				return true;
			}
		}
		return false;
	}

	/**
	 * @brief Retrieves a node with the given input from the graph.
	 * @param input A pointer to the input data to be checked.
	 * @return A shared pointer to the node if found, nullptr otherwise.
	 */
	std::shared_ptr<Node> getNode(double* input){
		for (int i = 0; i < nodes.size(); i++){
			if (nodes[i]->checkEqual(input)){
				return nodes[i];
			}
		}
		return nullptr;
	}

	/**
	 * @brief Prints the number of nodes in the graph.
	 */
	void printnumNodes(){
		printf("number of nodes in GRAPH: %zu\n", nodes.size());
	}

	/**
	 * @brief Prints the number of neighbors of a given node.
	 * @param node A shared pointer to the node whose neighbors are to be counted.
	 */
	void printnumNeighbors(std::shared_ptr<Node> node){
		if (edges.find(node) != edges.end()) {
			printf("number of neighbors: %zu\n", edges[node].size());
		}
		else{
			printf("Node has no neighbors found\n");
		}
	}

	/**
	 * @brief Checks if two nodes are neighbors.
	 * @param node1 A shared pointer to the first node.
	 * @param node2 A shared pointer to the second node.
	 * @return True if the nodes are neighbors, false otherwise.
	 */
	bool isNeighbor(std::shared_ptr<Node> node1, std::shared_ptr<Node> node2){
		if (edges.find(node1) != edges.end()) {
			for (int i = 0; i < edges[node1].size(); i++){
				if (edges[node1][i]->checkEqual(node2->angles)){
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * @brief Prints the neighbors of a given node.
	 * @param node A shared pointer to the node whose neighbors are to be printed.
	 */
	void printNeighbors(std::shared_ptr<Node> node){
		if (edges.find(node) != edges.end()) {
			for (int i = 0; i < edges[node].size(); i++){
				printf("neighbor %d: %f %f %f %f %f\n", i, edges[node][i]->angles[0], edges[node][i]->angles[1], edges[node][i]->angles[2], edges[node][i]->angles[3], edges[node][i]->angles[4]);
			}
		}
		else{
			printf("Node has no neighbors found\n");
		}
	}
};
struct Graph{
	std::vector<std::shared_ptr<Node>> nodes;
	std::unordered_map<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, NodePtrHash, NodePtrEqual> edges;

	Graph() {}

	void addNode(std::shared_ptr<Node> node) {
		nodes.push_back(node);
	}

	void addEdge(std::shared_ptr<Node> from, std::shared_ptr<Node> to) {
		//check if from and to are already in the graph
		
		if (edges.find(from) == edges.end()) {
			edges[from] = std::vector<std::shared_ptr<Node>>();
		}
		if (std::find(edges[from].begin(), edges[from].end(), to) == edges[from].end()) {
			edges[from].push_back(to);
		}
		//edges[from].push_back(to);

		if (edges.find(to) == edges.end()) {
			edges[to] = std::vector<std::shared_ptr<Node>>();
		}
		if (std::find(edges[to].begin(), edges[to].end(), from) == edges[to].end()) {
			edges[to].push_back(from);
		}
		//edges[to].push_back(from);
	}

	std::vector<std::shared_ptr<Node>> getNeighbors(std::shared_ptr<Node> node) {
		if (edges.find(node) != edges.end()) {
			return edges[node];
		}
		printf("no neighbors found\n");
		return std::vector<std::shared_ptr<Node>>();
	}

	bool inGraph(double* input){ 
		//return edges.find(node) != edges.end();
		for (int i = 0; i < nodes.size(); i++){
			if (nodes[i]->checkEqual(input)){
				return true;
			}
		}
		return false;
	}

	std::shared_ptr<Node> getNode(double* input){
		for (int i = 0; i < nodes.size(); i++){
			if (nodes[i]->checkEqual(input)){
				return nodes[i];
			}
		}
		return nullptr;
	}

	void printnumNodes(){
		printf("number of nodes in GRAPH: %zu\n", nodes.size());
	}

	void printnumNeighbors(std::shared_ptr<Node> node){
		if (edges.find(node) != edges.end()) {
			printf("number of neighbors: %zu\n", edges[node].size());
		}
		else{
			printf("Node has no neighbors found\n");
		}
	}

	bool isNeighbor(std::shared_ptr<Node> node1, std::shared_ptr<Node> node2){
		if (edges.find(node1) != edges.end()) {
			for (int i = 0; i < edges[node1].size(); i++){
				if (edges[node1][i]->checkEqual(node2->angles)){
					return true;
				}
			}
		}
		return false;
	}

	void printNeighbors(std::shared_ptr<Node> node){
		if (edges.find(node) != edges.end()) {
			for (int i = 0; i < edges[node].size(); i++){
				printf("neighbor %d: %f %f %f %f %f\n", i, edges[node][i]->angles[0], edges[node][i]->angles[1], edges[node][i]->angles[2], edges[node][i]->angles[3], edges[node][i]->angles[4]);
			}
		}
		else{
			printf("Node has no neighbors found\n");
		}
	}

};

/**
 * @brief Checks if a path exists between the start node and the goal node in the given graph.
 *
 * This function performs a breadth-first search (BFS) to determine if there is a path from the 
 * start node to the goal node. It uses a queue to explore nodes level by level and a set to 
 * keep track of visited nodes to avoid revisiting them.
 *
 * @param graph Pointer to the graph in which the search is performed.
 * @param startNode Shared pointer to the starting node.
 * @param goalNode Shared pointer to the goal node.
 * @return true if a path exists from the start node to the goal node, false otherwise.
 */
bool isPathExists(Graph* graph, std::shared_ptr<Node> startNode, std::shared_ptr<Node> goalNode) {
    if (startNode == nullptr || goalNode == nullptr) {
        return false;
    }

    std::queue<std::shared_ptr<Node>> toVisit;
    std::unordered_set<std::shared_ptr<Node>, NodePtrHash, NodePtrEqual> visited;

    toVisit.push(startNode);
    visited.insert(startNode);

    while (!toVisit.empty()) {
        std::shared_ptr<Node> currNode = toVisit.front();
        toVisit.pop();

        if (currNode->checkEqual(goalNode->angles)) {
			printf("Number of nodes visited: %ld\n", visited.size());
            return true;
        }

        std::vector<std::shared_ptr<Node>> neighbors = graph->getNeighbors(currNode);
        for (const auto& neighbor : neighbors) {
            if (visited.find(neighbor) == visited.end()) {
                toVisit.push(neighbor);
                visited.insert(neighbor);
            }
        }
    }
	printf("Number of nodes visited: %ld\n", visited.size());

    return false;
}

/**
 * @brief Builds a Probabilistic Roadmap (PRM) for motion planning.
 * 
 * This function constructs a PRM by sampling random configurations, adding them to the graph if they are valid,
 * and connecting them to their nearest neighbors within a specified range. It also ensures that the start and goal
 * nodes are connected to the graph.
 * 
 * @param graph Pointer to the Graph object where nodes and edges will be added.
 * @param tree Pointer to the Tree object used for nearest neighbor search.
 * @param map Pointer to the map data representing the environment.
 * @param x_size The width of the map.
 * @param y_size The height of the map.
 * @param num_samples The number of random samples to generate.
 * @param step_size The step size for sampling configurations.
 * @param dof The degrees of freedom of the robot.
 * @param range The range within which to search for nearest neighbors.
 * @param goal Pointer to the goal configuration.
 * @param start Pointer to the start configuration.
 */
void buildPRM(Graph* graph, Tree* tree, double* map, int x_size, int y_size, int num_samples, double step_size, int dof, double range, double* goal, double* start){ 
	
	int k=0;
	//printf("range %f\n", range);
	std::shared_ptr<Node> goalNode = graph->getNode(goal);
	std::shared_ptr<Node> startNode = graph->getNode(start);
	while(k<num_samples){
		k++;
		double* qrand = get_random_config(dof, map, x_size, y_size); //also check for valid configuration
		//printf("qrand angles %f %f %f %f %f\n", qrand[0], qrand[1], qrand[2], qrand[3], qrand[4]);
		if (graph->inGraph(qrand)){
			continue;
		}
		if (k%1000==0){
			printf("sample %d\n", k);
		}
		std::shared_ptr<Node> newNode = std::make_shared<Node>(nullptr, dof, qrand, nullptr, nullptr, nullptr, 0, 0, INT_MAX);
		graph->addNode(newNode);
		//printf("added node %d\n", k);
		tree->insert(newNode);
		vector<std::shared_ptr<Node>> neighbors;
		tree->nearestRangeNeighbor(qrand, tree->root, range, neighbors);
		if (neighbors.size() == 0){
			continue;
		}
		//printf("neighbors found %ld\n", neighbors.size());
		for (int i = 0; i < neighbors.size(); i++){
			if (neighbors[i]->checkValid(qrand, map, x_size, y_size, 7)){
				graph->addEdge(newNode, neighbors[i]);
				//newNode->cost = newNode->distance(goal);//////////////////////////
			}
		}
		if (k%1000==0){
		if(isPathExists(graph, startNode, goalNode)){
			printf("path exists\n");
			break;
		}
	}
		
	}
	vector<std::shared_ptr<Node>> neighbors;
	tree->nearestRangeNeighbor(goal, tree->root, range, neighbors);
		if (neighbors.size() == 0){
			printf("no neighbors found for goal\n");
			return;
		}
	
		for (int i = 0; i < neighbors.size(); i++){
			if (neighbors[i]->checkValid(goal, map, x_size, y_size)){
				
				graph->addEdge(goalNode, neighbors[i]);
				goalNode->cost = goalNode->distance(goal);
			}
		}
	
	neighbors.clear();
	tree->nearestRangeNeighbor(start, tree->root, range, neighbors);
	//printf("neighbors found for start %ld\n", neighbors.size());
	for (int i = 0; i < neighbors.size(); i++){
		if (neighbors[i]->checkValid(start, map, x_size, y_size)){
			//printf("angle %f %f %f %f %f\n", neighbors[i]->angles[0], neighbors[i]->angles[1], neighbors[i]->angles[2], neighbors[i]->angles[3], neighbors[i]->angles[4]);
			graph->addEdge(startNode, neighbors[i]);
			startNode->cost = startNode->distance(start);
		}
	}
	printf("PRM finished with %d samples\n", k);

	



}

bool isInPriorityQueue(const std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, CompareNode>& pq, const std::shared_ptr<Node>& node) {
    // Create a copy of the priority queue
    std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, CompareNode> pqCopy = pq;

    // Iterate through the elements of the copied priority queue
    while (!pqCopy.empty()) {
        if (pqCopy.top()->checkEqual(node->angles)) {
            return true;
        }
        pqCopy.pop();
    }
    return false;
}

void graphSearch(Graph* graph, double* start, double* goal, double*** plan, int* planlength) {
	// find the shortest path in the graph using priority queue

	std::shared_ptr<Node> startNode = graph->getNode(start);
	startNode->cost = 0;
	std::shared_ptr<Node> goalNode = graph->getNode(goal);
	goalNode->cost = INT_MAX;
	std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, CompareNode> open;
	std::unordered_set<std::shared_ptr<Node>, NodePtrHash, NodePtrEqual> closed;

	std::shared_ptr<Node> currNode= startNode ;
	open.push(currNode);
	int count=0;
	while (open.empty() == false){
		count++;
		
		currNode = open.top();
		open.pop();
		// if (currNode->checkEqual(goalNode->angles)){
		// 	printf("------------------goal found\n");
		// 	break;
		// }
		closed.insert(currNode);
		std::vector<std::shared_ptr<Node>> neighbors = graph->getNeighbors(currNode);
		if(currNode==nullptr){
			printf("curr node is null\n");
			break;
		}
		if(currNode->checkEqual(goalNode->angles)){
			printf("-----------------goal found\n");
			std::vector<std::shared_ptr<Node>> path;
			while (currNode != nullptr) {
				path.push_back(currNode);
				currNode = currNode->path_parent;
			}
        	std::reverse(path.begin(), path.end());

			*planlength = path.size();
			printf("plan length %d\n", *planlength);
			*plan = (double**)malloc(path.size() * sizeof(double*));
			for (int i = 0; i < path.size(); i++) {
				(*plan)[i] = (double*)malloc(goalNode->dof * sizeof(double));
				for (int j = 0; j < goalNode->dof; j++) {
					(*plan)[i][j] = path[i]->angles[j];
				}
			}
		
			break;
		}
		for (int i=0; i<neighbors.size(); i++){

			if (closed.find(neighbors[i]) != closed.end()){
				continue;
			}
			if (neighbors[i]->checkEqual(goalNode->angles)){
				//printf("found goal in neighbors\n");
			}
			double newCost = currNode->cost + currNode->distance(neighbors[i]->angles);
			if (newCost < neighbors[i]->cost){
				neighbors[i]->cost = newCost;
				neighbors[i]->path_parent = currNode;
				if(!isInPriorityQueue(open, neighbors[i])){
					open.push(neighbors[i]);
				}
				//open.push(neighbors[i]);
			}
		}


	}
	printf("Search finished\n");
	printf("count %d\n", count);
	printf("nodes in closed set %ld\n", closed.size());

}


bool verifyPath2P(double* p1, double* p2, double* map, int x_size, int y_size, int dof){
	double* new_angles = new double[dof];
	double* dist = new double[dof];
	double factor;
	for (int k=1; k<=5; k++){
		for (int i = 0; i < dof; i++){
			dist[i] = p2[i] - p1[i];
			factor = (double)k/(5+1);  // for k=2 factors: 0.33, 0.66, for k=3 factors: 0.25, 0.5, 0.75
			new_angles[i] = p1[i] + factor*dist[i];
		}
		if (!IsValidArmConfiguration(new_angles, dof, map, x_size, y_size)){
			delete[] new_angles;
			return false;
		}
	}
	delete[] new_angles;
	return true;
}

bool verifyPath(double*** plan, int length, double* map, int x_size, int y_size, int dof){
	//printf("initial node values %f %f %f %f %f\n", (*plan)[13][0], (*plan)[13][1], (*plan)[13][2], (*plan)[13][3], (*plan)[13][4]);
	
	//printf("Initial node valid: %d\n", IsValidArmConfiguration((*plan)[15], dof, map, x_size, y_size));
	for (int i = 1; i < length; i++){
		//printf("checking node %d\n", i);
		double* olnod = (*plan)[i-1];
		double* nod = (*plan)[i];
		if(nod== nullptr){
			printf("node is null\n");
			return false;
		}
		//printf("node i valid: %d\n", IsValidArmConfiguration(nod, dof, map, x_size, y_size));
		if (!IsValidArmConfiguration(nod, dof, map, x_size, y_size)){
			printf("invalid configuration at %d\n", i);
			return false;
		}
		if (!verifyPath2P(olnod, nod, map, x_size, y_size, dof)){
			printf("invalid path between %d and %d\n", i-1, i);
			return false;
		}
	}
	printf("valid path\n");
	return true;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                          DEFAULT PLANNER FUNCTION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//

void planner(
    double* map,
	int x_size,
	int y_size,
	double* armstart_anglesV_rad,
	double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
	double mind = INT_MAX;
    //for now just do straight interpolation between start and goal checking for the validity of samples
	srand(time(0));
	
	printf("5DOF planner\n");
	double* qrand;
	double* qrand2;
	for(int i=0; i<6; i++){
		qrand = get_random_config(5, map, x_size, y_size);
		qrand2 = get_random_config(5, map, x_size, y_size);
		printf(" %f %f %f %f %f\n", qrand[0], qrand[1], qrand[2], qrand[3], qrand[4]);
		printf(" %f %f %f %f %f\n", qrand2[0], qrand2[1], qrand2[2], qrand2[3], qrand2[4]);
		printf("\n");
	}
	printf("4DOF planner \n");
	for(int i=0; i<7; i++){
		qrand = get_random_config(4, map, x_size, y_size);
		qrand2 = get_random_config(4, map, x_size, y_size);
		printf(" %f %f %f %f\n", qrand[0], qrand[1], qrand[2], qrand[3]);
		printf(" %f %f %f %f\n", qrand2[0], qrand2[1], qrand2[2], qrand2[3]);
		printf("\n");
	}
	printf("3DOF planner \n");
	for(int i=0; i<13; i++){
		qrand = get_random_config(3, map, x_size, y_size);
		qrand2 = get_random_config(3, map, x_size, y_size);
		printf(" %f %f %f\n", qrand[0], qrand[1], qrand[2]);
		printf(" %f %f %f\n", qrand2[0], qrand2[1], qrand2[2]);
		printf("\n");
	}
	

    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return;
    }
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf) {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }    
    *planlength = numofsamples;
    
    return;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              RRT IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerRRT(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */
	srand(time(0));
	printf("start %f %f %f %f %f\n", armstart_anglesV_rad[0], armstart_anglesV_rad[1], armstart_anglesV_rad[2], armstart_anglesV_rad[3], armstart_anglesV_rad[4]);
	printf("Start valid %d\n", IsValidArmConfiguration(armstart_anglesV_rad, numofDOFs, map, x_size, y_size));
	printf("goal %f %f %f %f %f\n", armgoal_anglesV_rad[0], armgoal_anglesV_rad[1], armgoal_anglesV_rad[2], armgoal_anglesV_rad[3], armgoal_anglesV_rad[4]);
	std::shared_ptr<Node> start = std::make_shared<Node>(nullptr, numofDOFs, armstart_anglesV_rad);
	//printf("start node valid %d\n", IsValidArmConfiguration(start->angles, numofDOFs, map, x_size, y_size));
	Tree tree = Tree(start, numofDOFs, armstart_anglesV_rad);
	//printf("Tree root valid %d\n", IsValidArmConfiguration(tree.root->angles, numofDOFs, map, x_size, y_size));
	//printf("Tree initialized\n");
	tree.printTree();
	std::shared_ptr<Node> goal = std::make_shared<Node>(nullptr, numofDOFs, armgoal_anglesV_rad);
	//tree.insert(goal);
	//tree.extend(armgoal_anglesV_rad, 0.2, map, x_size, y_size);
	tree.build_RRT(armgoal_anglesV_rad, 10000000, 0.65, map, x_size, y_size);
    //planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
	bool found = tree.find(armgoal_anglesV_rad);
	printf("goal found %d\n", found);
	int depth;
	std::shared_ptr<Node> goalN;
	if (found){
		goalN = tree.findNode(armgoal_anglesV_rad);
		depth = goalN->pathlength; //tree.getDepth(tree.findNode(armgoal_anglesV_rad));
		printf("depth %d\n", depth);
	}
	else{
		double mind= INT_MAX;
		std::shared_ptr<Node> goalNN = tree.nearestNeighbor(armgoal_anglesV_rad, start, mind);
		printf("nearest neighbor %f %f %f %f %f\n", goalNN->angles[0], goalNN->angles[1], goalNN->angles[2], goalNN->angles[3], goalNN->angles[4]);
		printf("distance %f\n", goalNN->distance(armgoal_anglesV_rad));
		return;
	}
	
	//std::shared_ptr<Node> goalN = tree.findNode(armgoal_anglesV_rad);
	//printf("goalN angles %f %f %f %f %f\n", goalN->angles[0], goalN->angles[1], goalN->angles[2], goalN->angles[3], goalN->angles[4]);
	printf("depth %d\n", depth);
	tree.get_plan(goalN, plan, depth);
	//printf("plan index 0 valid %d\n", IsValidArmConfiguration((*plan)[0], numofDOFs, map, x_size, y_size));
	*planlength = depth+1;
	tree.printNumNodes();
	verifyPath(plan, *planlength, map, x_size, y_size, numofDOFs);
	double* goal_angles = (*plan)[0];
	//printf("goal angles %f %f %f %f %f\n", goal_angles[0], goal_angles[1], goal_angles[2], goal_angles[3], goal_angles[4]);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                         RRT CONNECT IMPLEMENTATION                                                //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerRRTConnect(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */
    srand(time(0));
	std::shared_ptr<Node> start = std::make_shared<Node>(nullptr, numofDOFs, armstart_anglesV_rad);
	std::shared_ptr<Node> goal = std::make_shared<Node>(nullptr, numofDOFs, armgoal_anglesV_rad);
	printf("start %f %f %f %f %f\n", armstart_anglesV_rad[0], armstart_anglesV_rad[1], armstart_anglesV_rad[2], armstart_anglesV_rad[3], armstart_anglesV_rad[4]);
	printf("goal %f %f %f %f %f\n", armgoal_anglesV_rad[0], armgoal_anglesV_rad[1], armgoal_anglesV_rad[2], armgoal_anglesV_rad[3], armgoal_anglesV_rad[4]);
	Tree tree = Tree(start, numofDOFs, armstart_anglesV_rad);
	Tree tree2 = Tree(goal, numofDOFs, armgoal_anglesV_rad);
	//printf("Tree initialized\n");
	double* link= RRTconnect(&tree, &tree2, armgoal_anglesV_rad, 900000, 0.55, map, x_size, y_size);
	if (link == nullptr){
		printf("No path found\n");
		return;
	}
	RRTconnectPath(&tree, &tree2, link, 0.4, map, x_size, y_size, plan, planlength);
	//printf("plan start %f %f %f %f %f\n", (*plan)[0][0], (*plan)[0][1], (*plan)[0][2], (*plan)[0][3], (*plan)[0][4]);
	//printf("plan end %f %f %f %f %f\n", (*plan)[*planlength-1][0], (*plan)[*planlength-1][1], (*plan)[*planlength-1][2], (*plan)[*planlength-1][3], (*plan)[*planlength-1][4]);
	printf("A plan found\n");
	printf("2 tree plan length %d\n", *planlength);
	tree.printNumNodes();
	tree2.printNumNodes();
	verifyPath(plan, *planlength, map, x_size, y_size, numofDOFs);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                           RRT STAR IMPLEMENTATION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerRRTStar(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */
	srand(time(0));
	std::shared_ptr<Node> start = std::make_shared<Node>(nullptr, numofDOFs, armstart_anglesV_rad);
	printf("start %f %f %f %f %f\n", armstart_anglesV_rad[0], armstart_anglesV_rad[1], armstart_anglesV_rad[2], armstart_anglesV_rad[3], armstart_anglesV_rad[4]);
	printf("goal %f %f %f %f %f\n", armgoal_anglesV_rad[0], armgoal_anglesV_rad[1], armgoal_anglesV_rad[2], armgoal_anglesV_rad[3], armgoal_anglesV_rad[4]);
	Tree tree = Tree(start, numofDOFs, armstart_anglesV_rad);
	printf("calling RRTSTAR\n");	
	tree.RRTstar(armgoal_anglesV_rad, 40000, 0.56, map, x_size, y_size);

	bool found = tree.find(armgoal_anglesV_rad);
	int depth;
	std::shared_ptr<Node> goalN;
	if (found){
		goalN = tree.findNode(armgoal_anglesV_rad);
		depth = goalN->pathlength; //tree.getDepth(tree.findNode(armgoal_anglesV_rad));
		printf("depth %d\n", depth);
	}
	else{
		double mind= INT_MAX;
		std::shared_ptr<Node> goalNN = tree.nearestNeighbor(armgoal_anglesV_rad, start, mind);
		printf("nearest neighbor %f %f %f %f %f\n", goalNN->angles[0], goalNN->angles[1], goalNN->angles[2], goalNN->angles[3], goalNN->angles[4]);
		printf("distance %f\n", goalNN->distance(armgoal_anglesV_rad));
		return;
	}
	
	//std::shared_ptr<Node> goalN = tree.findNode(armgoal_anglesV_rad);
	//printf("goalN angles %f %f %f %f %f\n", goalN->angles[0], goalN->angles[1], goalN->angles[2], goalN->angles[3], goalN->angles[4]);
	printf("depth %d\n", depth);
	tree.get_plan(goalN, plan, depth);
	*planlength = depth+1;
	printf("plan found\n");
	double* goal_angles = (*plan)[0];
	printf("goal angles %f %f %f %f %f\n", goal_angles[0], goal_angles[1], goal_angles[2], goal_angles[3], goal_angles[4]);
	verifyPath(plan, *planlength, map, x_size, y_size, numofDOFs);
	tree.printNumNodes();

}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              PRM IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerPRM(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */
	srand(time(0));
	std::shared_ptr<Node> start = std::make_shared<Node>(nullptr, numofDOFs, armstart_anglesV_rad, nullptr, nullptr, nullptr, 0, 0, 0);
	std::shared_ptr<Node> goal = std::make_shared<Node>(nullptr, numofDOFs, armgoal_anglesV_rad, nullptr, nullptr, nullptr, 0, 0, INT_MAX);
    printf("start %f %f %f %f %f\n", armstart_anglesV_rad[0], armstart_anglesV_rad[1], armstart_anglesV_rad[2], armstart_anglesV_rad[3], armstart_anglesV_rad[4]);
	printf("goal %f %f %f %f %f\n", armgoal_anglesV_rad[0], armgoal_anglesV_rad[1], armgoal_anglesV_rad[2], armgoal_anglesV_rad[3], armgoal_anglesV_rad[4]);
	Tree tree = Tree(start, numofDOFs, armstart_anglesV_rad);
	tree.insert(goal);
	double mind = INT_MAX;
	//Tree tree2 = Tree(goal, numofDOFs, armgoal_anglesV_rad);
	printf("\n");
	Graph graph;
	graph.addNode(start);
	graph.addNode(goal);
	//printf("dist start goal %f\n", start->distance(goal->angles));
	printf("\n");
	buildPRM(&graph, &tree, map, x_size, y_size, 12500, 0.5, numofDOFs, 0.75, armgoal_anglesV_rad, armstart_anglesV_rad);
	printf("------------------Neighbors built------------------\n");
	graph.printnumNeighbors(start);
	//graph.printNeighbors(start);
	//std::shared_ptr<Node> closestNeighbor = tree.nearestNeighbor2(armgoal_anglesV_rad, start, mind);
	//printf("closest neighbor %f %f %f %f %f\n", closestNeighbor->angles[0], closestNeighbor->angles[1], closestNeighbor->angles[2], closestNeighbor->angles[3], closestNeighbor->angles[4]);
	//printf("distance to goal %f\n", closestNeighbor->distance(armgoal_anglesV_rad));
	printf("----------nodes in graph \n");
	graph.printnumNodes();
	graphSearch(&graph, armstart_anglesV_rad, armgoal_anglesV_rad, plan, planlength);
	//printf("start and goal connected: %d\n", isPathExists(&graph, start, goal));	
	printf("goal node cost %f\n", goal->cost);
	double* goal_angles = (*plan)[0];
	printf("goal angles %f %f %f %f %f\n", goal_angles[0], goal_angles[1], goal_angles[2], goal_angles[3], goal_angles[4]);
	verifyPath(plan, *planlength, map, x_size, y_size, numofDOFs);
	tree.printNumNodes();
	//planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                MAIN FUNCTION                                                      //
//                                                                                                                   //
//*******************************************************************************************************************//

/** Your final solution will be graded by an grading script which will
 * send the default 6 arguments:
 *    map, numOfDOFs, commaSeparatedStartPos, commaSeparatedGoalPos, 
 *    whichPlanner, outputFilePath
 * An example run after compiling and getting the planner.out executable
 * >> ./planner.out map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 0 output.txt
 * See the hw handout for full information.
 * If you modify this for testing (e.g. to try out different hyper-parameters),
 * make sure it can run with the original 6 commands.
 * Programs that do not will automatically get a 0.
 * */
int main(int argc, char** argv) {
	double* map;
	int x_size, y_size;

    std::string mapDirPath = MAPS_DIR;
    std::string mapFilePath = mapDirPath + "/" + argv[1];
    std::cout << "Reading problem definition from: " << mapFilePath << std::endl;
	tie(map, x_size, y_size) = loadMap(mapFilePath);
	const int numOfDOFs = std::stoi(argv[2]);
	double* startPos = doubleArrayFromString(argv[3]);
	double* goalPos = doubleArrayFromString(argv[4]);
	int whichPlanner = std::stoi(argv[5]);

    std::string outputDir = OUTPUT_DIR;
	string outputFile = outputDir + "/" + argv[6];
	std::cout << "Writing solution to: " << outputFile << std::endl;

	if(!IsValidArmConfiguration(startPos, numOfDOFs, map, x_size, y_size)||
			!IsValidArmConfiguration(goalPos, numOfDOFs, map, x_size, y_size)) {
		throw runtime_error("Invalid start or goal configuration!\n");
	}

	///////////////////////////////////////
	//// Feel free to modify anything below. Be careful modifying anything above.

	double** plan = NULL;
	int planlength = 0;

    // Call the corresponding planner function
    if (whichPlanner == PRM)
    {
        plannerPRM(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else if (whichPlanner == RRT)
    {
		printf("RRT\n");
        plannerRRT(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else if (whichPlanner == RRTCONNECT)
    {
        plannerRRTConnect(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else if (whichPlanner == RRTSTAR)
    {
        plannerRRTStar(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else
    {	
        planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }

	//// Feel free to modify anything above.
	//// If you modify something below, please change it back afterwards as my 
	//// grading script will not work and you will recieve a 0.
	///////////////////////////////////////

    // Your solution's path should start with startPos and end with goalPos
    if (!equalDoubleArrays(plan[0], startPos, numOfDOFs) || 
    	!equalDoubleArrays(plan[planlength-1], goalPos, numOfDOFs)) {
		throw std::runtime_error("Start or goal position not matching");
	}

	/** Saves the solution to output file
	 * Do not modify the output log file output format as it is required for visualization
	 * and for grading.
	 */
	std::ofstream m_log_fstream;
	m_log_fstream.open(outputFile, std::ios::trunc); // Creates new or replaces existing file
	if (!m_log_fstream.is_open()) {
		throw std::runtime_error("Cannot open file");
	}
	m_log_fstream << mapFilePath << endl; // Write out map name first
	/// Then write out all the joint angles in the plan sequentially
	for (int i = 0; i < planlength; ++i) {
		for (int k = 0; k < numOfDOFs; ++k) {
			m_log_fstream << plan[i][k] << ",";
		}
		m_log_fstream << endl;
	}
}
