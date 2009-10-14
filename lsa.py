"""
Youen Peron
2009

LSA in python from Joseph Wilk implementation following Alex Thomo's tutorial.

http://blog.josephwilk.net/projects/latent-semantic-analysis-in-python.html

http://alexthomo.blogspot.com/2009/03/latent-semantic-analysis-tutorial.html

"""

from scipy import linalg,array,dot,mat,transpose
from math import *
from pprint import pprint
  
  
class LSA:
    """ Latent Semantic Analysis(LSA).
    Apply transforms to a document-term matrix to bring out latent relationships.
    These are found by analysing relationships between the documents and the terms they
    contain.
     """
  
  
    def __init__(self, matrix, termslabel, docslabel):
        self.matrix = array(matrix)
        self.termslabel = termslabel
        self.docslabel = docslabel
 
    def __repr__(self,):
        """ Make the matrix look pretty """
        stringRepresentation=""
  
        rows,cols = self.matrix.shape
      
        for row in xrange(0,rows):
            stringRepresentation += "["
  
            for col in xrange(0,cols):
                stringRepresentation+= "%+0.2f "%self.matrix[row][col]
            stringRepresentation += "]\n"
  
        return stringRepresentation
  
  
    
  
    def lsaTransform(self,k=2):
        """ Calculate SVD of objects matrix: U . SIGMA . VT = MATRIX
        Reduce the dimension of sigma by specified factor producing sigma'.
        Then dot product the matrices:  U . SIGMA' . VT = MATRIX'
        """
        rows,cols= self.matrix.shape

        if k < rows: #Its a valid reduction
  
            #Sigma comes out as a list rather than a matrix
            s,sigma,ut = linalg.svd(self.matrix)
      
            #Dimension reduction, build SIGMA'
            sigmak = linalg.diagsvd(sigma[:k],k,k)
            sk = transpose(transpose(s)[:k])
            utk = ut[:k] 
            #print utk
            #print sk.shape
            #print utk.shape

            return( dot(sk,sigmak), transpose(dot(sigmak,utk)))
        else:
         print "dimension reduction cannot be greater than %s" % rows

    def show(self):
        import matplotlib.pyplot as plt
        plt.axis([-4, 4, -4, 4])
        terms, docs = self.lsaTransform(2)
        for t, label in zip(terms,self.termslabel):
            plt.text(t[0],t[1],label)

        for d,label in  zip(docs,self.docslabel):
            plt.text(d[0],d[1],label)

        plt.show()

 
  
if __name__ == '__main__':
  
    #Example document-term matrix
    # Vector  d1,  d2,  d3,  d4,  d5
    matrix = [[1.0, 0.0, 1.0, 0.0, 0.0], # romeo
             [1.0, 1.0, 0.0, 0.0, 0.0], # juliet
             [0.0, 1.0, 0.0, 0.0, 0.0], # happy
             [0.0, 1.0, 1.0, 0.0, 0.0], # dagger
             [0.0, 0.0, 0.0, 1.0, 0.0], # live
             [0.0, 0.0, 1.0, 1.0, 0.0], # die
             [0.0, 0.0, 0.0, 1.0, 0.0], # free
             [0.0, 0.0, 0.0, 1.0, 1.0]] # new-hampshire
 
    #Create
    termslabel = ["romeo","juliet","happy","dagger","live","die","free","new-hampshire"]
    docslabel = ["d1","d2","d3","d4","d5"]
    lsa = LSA(matrix,termslabel, docslabel )
    
    terms,docs = lsa.lsaTransform(2)
    print "Terms"
    print terms
    print "Docs"
    print docs
    
    lsa.show()

