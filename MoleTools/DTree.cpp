//Sean M. Law

#include "DTree.hpp"

DTree::DTree(){
	root=NULL;
}

DTree::~DTree(){
	delDTree();
}

void DTree::delDTree(DTreeNode *leaf){
	if (leaf != NULL){
		delDTree(leaf->left);
		delDTree(leaf->right);
		delete leaf;
	}
}

/*
void DTree::addDTree(double key, DTreeNode *leaf){
	if (key < leaf ->key_value){
		if (leaf->left != NULL){
			addDTree(key, leaf->left);
		}
		else{
			leaf->left=new DTreeNode;
			leaf->left->key_value=key;
			leaf->left->left=NULL; //Sets left child of child node to NULL
			leaf->left->right=NULL; //Sets right child of child node to NULL
		}
	}
	else{
		if (leaf->right != NULL){
			addDTree(key, leaf->right);
		}
		else{
			leaf->right=new DTreeNode;
			leaf->right->key_value=key;
			leaf->right->left=NULL; //Sets left child of child node to NULL
			leaf->right->right=NULL; //Sets right child of child node to NULL
		}
	}
}

DTreeNode* DTree::searchDTree(double key, DTreeNode *leaf){
	if (leaf != NULL){
		if (key == leaf->key_value){
			return leaf;
		}
		else if (key < leaf->key_value){
			return searchDTree(key, leaf->left);
		}
		else{
			return searchDTree(key, leaf->right);
		}
	}
	else{
		return NULL;
	}
}
*/

void DTree::addDTree(double key, DTreeNode *leaf, unsigned int index, std::string classin){
  if (leaf->left != NULL){
    addDTree(key, leaf->left, index, classin);
  }
  else{
    leaf->left=new DTreeNode;
    leaf->left->key_value=key;
		leaf->left->inx=index;
    leaf->left->left=NULL; //Sets left child of child node to NULL
    leaf->left->right=NULL; //Sets right child of child node to NULL
		leaf->left->cls=classin; //The class of this leaf node
  }
  if (leaf->right != NULL){
    addDTree(key, leaf->right, index, classin);
  }
  else{
    leaf->right=new DTreeNode;
    leaf->right->key_value=key;
		leaf->right->inx=index;
    leaf->right->left=NULL; //Sets left child of child node to NULL
    leaf->right->right=NULL; //Sets right child of child node to NULL
		leaf->right->cls=classin; //The class of this leaf node
  }
}

std::string DTree::getDTreeClass(DTreeNode *leaf, std::vector<std::vector<double> > fin){
	if (leaf->left != NULL){
		return getDTreeClass(leaf->left, fin);
	}
	else{
		return leaf->cls;
	}
}

/*
void DTree::addDTree(double key){
	//Public version of the addDTree function
	//Takes care of when the root node has not been initialized
	if (root != NULL){
		addDTree(key, root);
	}
	else{
		root=new DTreeNode;
		root->key_value=key;
		root->left=NULL;
		root->right=NULL;
	}
}

DTreeNode* DTree::searchDTree(double key){
	//Public version of the searchDTree function
	//Starts search recursion starting at root node
	return searchDTree(key, root);
}
*/

void DTree::addDTree(double key, unsigned int index, std::string classin){
  //Public version of the addDTree function
  //Takes care of when the root node has not been initialized
  if (root != NULL){
    addDTree(key, root, index, classin);
  }
  else{
    root=new DTreeNode;
    root->key_value=key;
		root->inx=index;
    root->left=NULL;
    root->right=NULL;
		root->cls=classin;
  }
}



void DTree::delDTree(){
	//Public version of the delDTree function
	//Starts delete recursion starting at root node
  delDTree(root);
}

DTreeNode* DTree::getDTreeRoot(){
	if (root != NULL){
		return root;
	}
	else{
		return NULL;
	}
}

std::string DTree::getDTreeClass(std::vector<std::vector<double> > fin){
	//Public version of the getDTreeClass function
	if (root != NULL){
		return getDTreeClass(root, fin);
	}
	else{
		return "?";	
	}
}

DTreeNode* DTreeNode::getDTreeNodeLeft(){
	return this->left;
}

DTreeNode* DTreeNode::getDTreeNodeRight(){
	return this->right;
}
