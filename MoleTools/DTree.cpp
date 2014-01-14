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

void DTree::genDTree(DTreeNode *&leaf, std::vector<std::string> &t, unsigned int &inx, std::string delim){
	//Pre-order string expected
	//Example: 1.0:1 2.0:2 3.0:3 4.0:4 H C 5.0:5 E H 6.0:6 C E 7.0:7 8.0:8 9.0:9 E C 10.0:10 H E 11.0:11 H E 
	//Value:Index or Class

  std::vector<std::string> s;
  unsigned int i;
  double d;

  if (delim.compare(0,1," " ) == 0){
    std::cerr << "Error: Delimiter cannot be whitespace!" << std::endl;
    return;
  }

  if (inx < t.size()){
    leaf=new DTreeNode;
    leaf->left=NULL;
    leaf->right=NULL;

    Misc::splitStr(Misc::trim(t.at(inx)), delim, s, false);
		if (s.size() == 1){
      //std::cerr << t.at(inx) << std::endl;
      leaf->cls=s.at(0);
      return;
    }
    else if (s.size() == 2){
      //std::cerr << t.at(inx) << std::endl;
      std::stringstream(s.at(0)) >> d;
      leaf->key_value=d;
      std::stringstream(s.at(1)) >> i;
      leaf->inx=i-1;
      this->genDTree(leaf->left, t, ++inx, delim);
			this->genDTree(leaf->right, t, ++inx, delim);
    }
    else{
      //I shouldn't be here
    }
  }
  else{
    //Do nothing
  }
}

std::string DTree::getDTreeClass(DTreeNode *leaf, const std::vector<double> &fin){
	if (leaf->left == NULL && leaf->right == NULL){
		return leaf->cls;
	}
	else if (leaf->left == NULL || leaf->right == NULL){
		if (leaf->left != NULL){
			return getDTreeClass(leaf->left, fin);
		}
		else{
			return getDTreeClass(leaf->right, fin);
		}
	}
	else{
		if (leaf->inx >= fin.size()){
			std::cerr << "Error: Missing feature (inx = " << leaf->inx << " )" << std::endl;
			return "";
		}
		else if (fin.at(leaf->inx) < leaf->key_value){
			return getDTreeClass(leaf->left, fin);
		}
		else{
			return getDTreeClass(leaf->right, fin);
		}
	}
}

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

void DTree::genDTree(std::vector<std::string> &t, std::string delim){
	//Public version of the genDTree function
	//Pre-order string expected
	//Example: 1.0:1 2.0:2 3.0:3 4.0:4 H C 5.0:5 E H 6.0:6 C E 7.0:7 8.0:8 9.0:9 E C 10.0:10 H E 11.0:11 H E
	//Value:Index or Class

	unsigned int inx;
	std::vector<std::string> s;
	unsigned int i; 
	double d;

	inx=0;

	if (delim.compare(0,1," " ) == 0){
		std::cerr << "Error: Delimiter cannot be whitespace!" << std::endl;
		return;
	}

	if (inx < t.size()){
		root=new DTreeNode;
    root->left=NULL;
    root->right=NULL;
	
		Misc::splitStr(Misc::trim(t.at(inx)), delim, s, false);
		if (s.size() == 1){
			//std::cerr << t.at(inx) << std::endl;
			root->cls=s.at(0);
			return;
		}
		else if (s.size() == 2){
			//std::cerr << t.at(inx) << std::endl;
			std::stringstream(s.at(0)) >> d;
			root->key_value=d;
			std::stringstream(s.at(1)) >> i;
			root->inx=i-1;
			this->genDTree(root->left, t, ++inx, delim);
			this->genDTree(root->right, t, ++inx, delim);
		}
		else{
			//I shouldn't be here
		}
	}
	else{
		//Do nothing
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

std::string DTree::getDTreeClass(const std::vector<double> &fin){
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
