#install packages
import synapseclient 
import synapseutils 
import os

#access UKBB protein data from synapse
syn = synapseclient.Synapse() 
token = "" ##enter personal token from synapse profile
syn.login(authToken=token)
destination_folder = "ukbb_pqtl"
files = synapseutils.syncFromSynapse(syn, 'syn51365303', path=destination_folder, ifcollision='keep.local')
