#!/usr/bin/env python
# encoding: utf-8
"""
SimpleHMM.py
A basic implementation of hidden Markov modelling.  Includes Baum-Welch EM to estimate transition matrices

Created by Joshua Shapiro on 2008-11-11.
"""

import sys
import os
from numpy import *



class HMM(object):
  """Object that implements a hidden markov model, given a matrix of likelihoods."""
  def __init__(self, likeMatrix, transMatrix):
    super(HMM, self).__init__()
    self.like = array(likeMatrix, copy = False)
    self._trans = matrix(transMatrix)
    self.nstates = self.like.shape[0]
    self.nsites = self.like.shape[1]
    # check that matrices are the correct dimensions
    if (self._trans.shape[0] != self.nstates 
      or self._trans.shape[1] != self.nstates):
       raise ValueError, "Number of states and number of rows in likelihood matrix must agree"
    # check that transition matrix is well formed( columns add to 1)
    for col in self._trans.T:
      if sum(col) - 1.0 > 1e-7:
        raise ValueError, "The sum of the transition matrix columnsis not 1"
    self._forward = None
    self._backward = None
    self._scale = zeros(self.nsites)
    self._totalLike = None
    self._posterior = None
  
  def getTrans(self):
    return self._trans
  
  def setTrans(self, transMatrix):
    """Sets a new Transition matrix and resets all calculations"""
    self._trans = matrix(transMatrix)
    self._forward = None
    self._backward = None
    self._scale = zeros(self.nsites)
    self._totalLike = None
    self._posterior = None
  
  def getScale(self):
    if self._scale[0] == 0.:
      self.calcForward()
    return self._scale
  
  def getTotal(self):
    if self._totalLike == None:
      self.calcForward()
    return self._totalLike
  
  
  def calcForward(self):
     """Forward HMM matrix likelihoods"""
     if self._forward == None:
       fLikes =  zeros((self.nstates, self.nsites))
       #first position with equal initial probs
       fLikes[:,0] = (self.like[:,0]) /self.nstates
       self._scale[0] = sum(fLikes[:,0])
       fLikes[:,0] /= self._scale[0]
       self._totalLike = log(self._scale[0])
       for i in range(self.nsites - 1):
         fLikes[:,i+1] = array(fLikes[:,i] * self._trans.T) * self.like[:,i+1]
         self._scale[i+1] = sum(fLikes[:,i+1])
         fLikes[:,i+1] /= self._scale[i+1]
         self._totalLike += log(self._scale[i+1])
       self._forward = fLikes
     return self._forward

  def calcBackward(self):
    """Calculate Backward HMM matrix likelihoods"""
    if self._backward == None:
      scale = self.scale
      bLikes = zeros((self.nstates, self.nsites))
      #last position
      bLikes[:,-1] = 1.0/scale[-1]
      for i in reversed(range(0, self.nsites-1)):
        bLikes[:,i] = array((bLikes[:,i+ 1] * self.like[:,i+ 1])  * self._trans) / scale[i]
      self._backward = bLikes
    return self._backward

  def calcPosterior(self):
    """Calculate the posterior probabilities for each trait"""
    if self._posterior == None:
      self._posterior = (self.forward * self.backward) * self.scale
    return self._posterior

  def baumWelch(self, nIter = 100, verbose = False):
    """Use Baum-Welch EM algorithm to update transition matrix.  Returns log likelihood of data for each iteration"""
    if verbose:
      print "Iteration\tLog Likelihood"
      print "0\t", self.logLike
    progress = [self.logLike]
    for i in range(nIter):
      self.transition = (dot(self.backward[:,1:] * self.like[:,1:], self.forward[:,:(self.nsites-1)].T) * array(self._trans) 
                     / sum(self.posterior[:,:self.nsites-1], axis = 1) )
      if verbose:
        print i+1, "\t" , self.logLike
      progress.append(self.logLike)
    return progress

  def calcViterbi(self):
    """Calculate the Viterbi decoding"""
    path = [[] for i in range(self.nstates)]
    lastProb = [log(1.0) for i in range(self.nstates)]
   
    for i in range(self.nsites):
      tempPath = path
      tempProb = lastProb
      for state in range(self.nstates):
        transProbs = [lastProb[last] + log(self._trans[state, last]) for last in range(self.nstates)]
        tempProb[state] = max(transProbs) + log(self.like[state, i])
        tempPath[state] = path[transProbs.index(max(transProbs))] + [state]
      path = tempPath
      lastProb = tempProb
    maxState = lastProb.index(max(lastProb))
    return (path[maxState])
      
                     
  
  #properties
  scale = property(fget = getScale, doc = "Scaling parameter for likelihoods")
  logLike = property(fget = getTotal,  doc ="Total log likelihood of data")
  forward = property(fget = calcForward,  doc ="Forward likelihood matrix")
  backward = property(fget = calcBackward,  doc ="Backward likelihood matrix")
  posterior = property(fget = calcPosterior,  doc ="Posterior probabilities for the states of the HMM")
  transition = property(fget = getTrans, fset = setTrans, doc = "Transition matrix")  
  viterbi = property(fget = calcViterbi, doc = "Viterbi decoding")
    



if __name__ == '__main__':
  #likes = array([[.6,.6,.2,.01,.2,.2,.2,.9,.6,.6],
  #              [.1,.1,.1,.9,.6,.6,.6,.01,.1,.1]])
  likes = array([[.9,.9,.9,.9,.9,.9,.9,.9,.9,.9],
                [.1,.1,.1,.1,.1,.1,.1,.1,.1,.1]])
  
  
  trans= matrix([[0.9,0.1],
                 [0.1,0.9]])
  test = HMM(likes, trans)
  print test.forward
  print test.backward
  print test.scale
  print ""
  print test.posterior
  print "-----------------------------"
  print test.viterbi
  print test.transition
  test.baumWelch(10, True)
  #print test.posterior
  #print test.transition
  #print test.viterbi
  