//Sean M. Law

#include "BTree.hpp"

#include <cstdlib>

BTree::BTree(){
	root=NULL;
}

BTree::~BTree(){
	delBTree();
}

void BTree::delBTree(BTreeNode *node){
	if (node != NULL){
		delBTree(node->left);
		delBTree(node->right);
		delete node;
	}
}

void BTree::addBTree(int key, BTreeNode *node){
	if (key < node ->key_value){
		if (node->left != NULL){
			addBTree(key, node->left);
		}
		else{
			node->left=new BTreeNode;
			node->left->key_value=key;
			node->left->left=NULL; //Sets left child of child node to NULL
			node->left->right=NULL; //Sets right child of child node to NULL
		}
	}
	else{
		if (node->right != NULL){
			addBTree(key, node->right);
		}
		else{
			node->right=new BTreeNode;
			node->right->key_value=key;
			node->right->left=NULL; //Sets left child of child node to NULL
			node->right->right=NULL; //Sets right child of child node to NULL
		}
	}
}

BTreeNode* BTree::searchBTree(int key, BTreeNode *node){
	if (node != NULL){
		if (key == node->key_value){
			return node;
		}
		else if (key < node->key_value){
			return searchBTree(key, node->left);
		}
		else{
			return searchBTree(key, node->right);
		}
	}
	else{
		return NULL;
	}
}

void BTree::addBTree(int key){
	//Public version of the addBTree function
	//Takes care of when the root node has not been initialized
	if (root != NULL){
		addBTree(key, root);
	}
	else{
		root=new BTreeNode;
		root->key_value=key;
		root->left=NULL;
		root->right=NULL;
	}
}

BTreeNode* BTree::searchBTree(int key){
	//Public version of the searchBTree function
	//Starts search recursion starting at root node
	return searchBTree(key, root);
}

void BTree::delBTree(){
	//Public version of the delBTree function
	//Starts delete recursion starting at root node
  delBTree(root);
}

