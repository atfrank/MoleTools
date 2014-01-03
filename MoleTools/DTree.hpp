//Sean M. Law

#ifndef DTREE_H
#define DTREE_H

#include <cstdlib>

struct DTreeNode{
  int key_value;
	DTreeNode *left;
	DTreeNode *right;
};


class DTree {
	private:
		void delDTree(DTreeNode *leaf);
		void addDTree(int key, DTreeNode *leaf);
		DTreeNode *searchDTree(int key, DTreeNode *leaf);
		DTreeNode *root;

	public:
		DTree();
		~DTree();

		void delDTree();
		void addDTree(int key);
		DTreeNode *searchDTree(int key);
};

#endif
