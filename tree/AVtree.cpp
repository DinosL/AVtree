#include "../def.h"

int Node::max(int a, int b)
{
    return (a > b)? a : b;
}

Node* Node::newNode(float epsilon, int start, int end, int id, Node *prevNode, int leftchild)
{
    Node* node = new Node();
    node->epsilon = epsilon;
    node->rec = NULL;
    node->id = id;
    node->start = start;
    node->end = end;
    node->left = NULL;
    node->right = NULL;
    node->leftchild = leftchild;
    node->prev = prevNode;
    return(node);
}

Node* Node::newNode(float epsilon, int start, int end, int id, int leftchild)
{
    Node* node = new Node();
    node->epsilon = epsilon;
    node->rec = NULL;
    node->id = id;
    node->start = start;
    node->end = end;
    node->left = NULL;
    node->right = NULL;
    node->leftchild = leftchild;
    return(node);
}


Node* Node::newNode(float epsilon, int start, int end, int id, Node *prevNode)
{
    Node* node = new Node();
    node->epsilon = epsilon;
    node->rec = NULL;
    node->id = id;
    node->start = start;
    node->end = end;
    node->left = NULL;
    node->right = NULL;
    node->prev = prevNode;
    return(node);
}

Node* Node::newNode(float epsilon, int start, int end, int id)
{
    Node* node = new Node();
    node->epsilon = epsilon;
    node->rec = NULL;
    node->id = id;
    node->start = start;
    node->end = end;
    node->left = NULL;
    node->right = NULL;
    return(node);
}

Node* Node::insert(Node* node, float epsilon, int start, int end, int id, Node* prevNode, int leftchild)
{
    if (node == NULL)
    {
        return(newNode(epsilon, start, end, id, prevNode,leftchild));
    }

    return node;
}


Node* Node::insert(Node* node, float epsilon, int start, int end, int id, Node* prevNode)
{
    if (node == NULL)
    {
        return(newNode(epsilon, start, end, id, prevNode));
    }
    return node;
}


Node* Node::insert(Node* node, float epsilon, int start, int end, int id, int leftchild)
{
    if (node == NULL)
    {
        return(newNode(epsilon, start, end, id, leftchild));
    }

    return node;
}

Node* Node::insert(Node* node, float epsilon, int start, int end, int id)
{
    if (node == NULL)
    {
        return(newNode(epsilon, start, end, id));
    }

    return node;
}

Node* Node::insertLeft(Node* node, float epsilon, int start, int end, int id, int leftchild)
{
    if (node == NULL)
    {
        return(newNode(epsilon, start, end, id, node, leftchild));
    }

    node->left = insert(node->left, epsilon, start, end, id, node, leftchild);
    return node;
}

Node* Node::insertRight(Node* node, float epsilon, int start, int end, int id, int leftchild)
{
    if (node == NULL)
    {
        return(newNode(epsilon, start, end, id, node, leftchild));
    }

    node->right = insert(node->right, epsilon, start, end, id, node, leftchild);
    return node;
}

Node* Node::insertLeft(Node* node, float epsilon, int start, int end, int id)
{
    if (node == NULL)
    {
        return(newNode(epsilon, start, end, id, node));
    }

    node->left = insert(node->left, epsilon, start, end, id, node);
    return node;
}

Node* Node::insertRight(Node* node, float epsilon, int start, int end, int id)
{
    if (node == NULL)
    {
        return(newNode(epsilon, start, end, id, node));
    }

    node->right = insert(node->right, epsilon, start, end, id, node);
    return node;
}


void Node::preOrder(Node *root)
{
    if(root != NULL)
    {
        cout << "[" << root->id << "," << root->start << "," << root->end << "] - " << root->epsilon << endl ;
        preOrder(root->left);
        preOrder(root->right);
    }
}

