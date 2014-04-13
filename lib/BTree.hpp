//Sean M. Law

#ifndef BTREE_H
#define BTREE_H

struct BTreeNode{
  int key_value;
	BTreeNode *left;
	BTreeNode *right;
};


class BTree {
	private:
		void delBTree(BTreeNode *node);
		void addBTree(int key, BTreeNode *node);
		BTreeNode *searchBTree(int key, BTreeNode *node);
		BTreeNode *root;

	public:
		BTree();
		~BTree();

		void delBTree();
		void addBTree(int key);
		BTreeNode *searchBTree(int key);
};

#endif
