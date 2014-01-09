//Sean M. Law

#ifndef DTREE_H
#define DTREE_H

#include "Misc.hpp"
#include <cstdlib>
#include <string>
#include <vector>

struct DTreeNode{
  double key_value;
	unsigned int inx; //1-D vector<double> index
	DTreeNode *left;
	DTreeNode *right;
	std::string cls;

	DTreeNode* getDTreeNodeLeft();
	DTreeNode* getDTreeNodeRight();

};


class DTree {
	private:

		void delDTree(DTreeNode *leaf);
/*
		void addDTree(double key, DTreeNode *leaf);
		DTreeNode* searchDTree(double key, DTreeNode *leaf);
*/
		DTreeNode* root;
		
		void addDTree(double key, DTreeNode *leaf, unsigned int index, std::string classin="");
		std::string getDTreeClass(DTreeNode *leaf, std::vector<std::vector<double> > fin);

	public:
		DTree();
		~DTree();

		void delDTree();
/*
		void addDTree(double key);
		DTreeNode* searchDTree(double key);
*/

		void addDTree(double key, unsigned int index, std::string classin="");
		DTreeNode* getDTreeRoot();
		std::string getDTreeClass(std::vector<std::vector<double> > fin);
};

#endif


