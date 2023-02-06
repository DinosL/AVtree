#include "../def.h"

Node* Node::newNodeString(float epsilon, int start, int end, int id, Node *prevNode, int leftchild)
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

Node* Node::newNodeString(float epsilon, int start, int end, int id, int leftchild)
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


Node* Node::newNodeString(float epsilon, int start, int end, int id, Node *prevNode)
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

Node* Node::newNodeString(float epsilon, int start, int end, int id)
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

Node* Node::insertString(Node* node, float epsilon, int start, int end, int id, Node* prevNode, int leftchild)
{
    if (node == NULL)
    {
        return(newNodeString(epsilon, start, end, id, prevNode,leftchild));
    }

    return node;
}


Node* Node::insertString(Node* node, float epsilon, int start, int end, int id, Node* prevNode)
{
    if (node == NULL)
    {
        return(newNodeString(epsilon, start, end, id, prevNode));
    }

    return node;
}


Node* Node::insertString(Node* node, float epsilon, int start, int end, int id, int leftchild)
{
    if (node == NULL)
    {
        return(newNodeString(epsilon, start, end, id, leftchild));
    }

    return node;
}

Node* Node::insertString(Node* node, float epsilon, int start, int end, int id)
{
    if (node == NULL)
    {
        return(newNodeString(epsilon, start, end, id));
    }

    return node;
}

Node* Node::insertLeftString(Node* node, float epsilon, int start, int end, int id, int leftchild)
{
    if (node == NULL)
    {
        return(newNodeString(epsilon, start, end, id, node, leftchild));
    }

    node->left = insertString(node->left, epsilon, start, end, id, node, leftchild);
    return node;
}

Node* Node::insertRightString(Node* node, float epsilon, int start, int end, int id, int leftchild)
{
    if (node == NULL)
    {
        return(newNodeString(epsilon, start, end, id, node, leftchild));
    }

    node->right = insertString(node->right, epsilon, start, end, id, node, leftchild);
    return node;
}

Node* Node::insertLeftString(Node* node, float epsilon, int start, int end, int id)
{
    if (node == NULL)
    {
        return(newNodeString(epsilon, start, end, id, node));
    }

    node->left = insertString(node->left, epsilon, start, end, id, node);
    return node;
}

Node* Node::insertRightString(Node* node, float epsilon, int start, int end, int id)
{
    if (node == NULL)
    {
        return(newNodeString(epsilon, start, end, id, node));
    }

    node->right = insertString(node->right, epsilon, start, end, id, node);
    return node;
}


void Node::preOrderString(Node *root)
{
    if(root != NULL)
    {
        cout << "[" << root->id << "," << root->start << "," << root->end << "] - " << root->epsilon << endl ;
        preOrderString(root->left);
        preOrderString(root->right);
    }
}

