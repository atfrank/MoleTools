//Sean M. Law

#ifndef DTREE_H
#define DTREE_H

#include <cstdlib>
#include <string>

struct DTreeNode{
  double key_value;
	unsigned int col;
	DTreeNode *left;
	DTreeNode *right;
	std::string SSE;
};


class DTree {
	private:
		void delDTree(DTreeNode *leaf);
		void addDTree(double key, DTreeNode *leaf);
		DTreeNode* searchDTree(double key, DTreeNode *leaf);
		DTreeNode* root;

	public:
		DTree();
		~DTree();

		void delDTree();
		void addDTree(double key);
		DTreeNode* searchDTree(double key);
};

#endif
