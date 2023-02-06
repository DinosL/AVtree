
class Node
{
    public:
	    int value;
	    int start, end;
		int id;
	    Node *left;
	    Node *right;
		Node *prev;
	    char *rec;
	    float epsilon;
		int leftchild;


	    Node* newNodeString(float epsilon, int start, int end, int id, int leftchild);
		Node* newNodeString(float epsilon, int start, int end, int id, Node *prev, int leftchild);
		Node* newNodeString(float epsilon, int start, int end, int id);
		Node* newNodeString(float epsilon, int start, int end, int id, Node *prev);
	    Node* insertString(Node* node, float epsilon, int start, int end, int id, int leftchild);
		Node* insertString(Node* node, float epsilon, int start, int end, int id, Node *prevNode, int leftchild);
	    Node* insertRightString(Node* node, float epsilon, int start, int end, int id, int leftchild);
	    Node* insertLeftString(Node* node, float epsilon, int start, int end, int id, int leftchild);
		Node* insertString(Node* node, float epsilon, int start, int end, int id);
		Node* insertString(Node* node, float epsilon, int start, int end, int id, Node *prevNode);
	    Node* insertRightString(Node* node, float epsilon, int start, int end, int id);
	    Node* insertLeftString(Node* node, float epsilon, int start, int end, int id);
	    void preOrderString(Node *root);
};