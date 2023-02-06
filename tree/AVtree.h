class Node
{
    public:
	    int value;
	    int start, end;
		int id;
	    Node *left;
	    Node *right;
		Node *prev;
	    float *rec;
	    float epsilon;
		int leftchild;


	    int max(int a, int b);
	    Node* newNode(float epsilon, int start, int end, int id, int leftchild);
		Node* newNode(float epsilon, int start, int end, int id, Node *prev, int leftchild);
		Node* newNode(float epsilon, int start, int end, int id);
		Node* newNode(float epsilon, int start, int end, int id, Node *prev);
	    Node* insert(Node* node, float epsilon, int start, int end, int id, int leftchild);
		Node* insert(Node* node, float epsilon, int start, int end, int id, Node *prevNode, int leftchild);
	    Node* insertRight(Node* node, float epsilon, int start, int end, int id, int leftchild);
	    Node* insertLeft(Node* node, float epsilon, int start, int end, int id, int leftchild);

		Node* insert(Node* node, float epsilon, int start, int end, int id);
		Node* insert(Node* node, float epsilon, int start, int end, int id, Node *prevNode);
	    Node* insertRight(Node* node, float epsilon, int start, int end, int id);
	    Node* insertLeft(Node* node, float epsilon, int start, int end, int id);
	    void preOrder(Node *root);
};