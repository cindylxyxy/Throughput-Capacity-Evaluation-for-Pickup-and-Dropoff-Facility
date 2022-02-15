########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file defines useful objects.
### WARNING ###
# The following section of the file should not be changed even if you know exactly what you are trying to do.
# WARNING: Unintended changes to this section can have unexpected effects to the numerical results.

########################################################################################################################################################
def sortHelper(headNode, tlen1, tlen2):

	curr = headNode
	print ('starting...')
	while curr.next:
		print (curr.data)
		curr = curr.next
	print (curr.data)
	print ('target length is:', tlen1, tlen2 )

	if headNode != None and headNode.next != None:
		prev = headNode
		curr = headNode.next
		assert tlen1 >= 1
		assert tlen2 >= 1
		assert tlen1 >= tlen2
		assert tlen1 - tlen2 <= 1
		for _ in range(tlen1 - 1):
			prev = curr
			curr = curr.next
		prev.next = None
			
		curr1 = sortHelper(headNode, tlen1 - tlen1//2, tlen1//2)
		curr2 = sortHelper(curr, tlen2 - tlen2//2, tlen2//2)

		if curr1.data <= curr2.data:
			headNode = curr1
			curr1 = curr1.next
		else:
			headNode = curr2
			curr2 = curr2.next
		curr = headNode 

		while curr1 is not None and curr2 is not None:
			if curr1.data <= curr2.data:
				curr.next = curr1
				curr1 = curr1.next
			else:
				curr.next = curr2
				curr2 = curr2.next
			curr = curr.next

		if curr1 is not None:
			assert curr2 is None
			curr.next = curr1
		if curr2 is not None:
			assert curr1 is None
			curr.next = curr2
	return headNode


class Node:

	def __init__(self, data = None):
		self.data = data
		self.next = None

class SLinkedList:

	def __init__(self):
		self.head = None
		self.length = 0
		
	def __len__(self):
		return self.length

	def print(self):
		curr = self.head
		while (curr):
			print (curr.data)
			curr = curr.next

	def isEmpty(self):
		assert (self.head is None) == (self.length == 0)
		return (self.head is None)

	def addHead(self, new):
		NewNode = Node(new)
		NewNode.next = self.head
		self.head = NewNode
		self.length += 1

	def addEnd(self, new):
		NewNode = Node(new)
		if self.head is None:
			self.head = NewNode
			return
		curr = self.head
		while (curr.next):
			curr = curr.next
		curr.next = NewNode
		self.length += 1

	def insert(self, loc, new):
		if loc == 0:
			self.addHead(new)
		elif loc >= self.length:
			self.addEnd(new)
		else:
			NewNode = Node(new)
			prev = self.head
			for _ in range(loc - 1):
				prev = prev.next
			curr = prev.next
			prev.next = NewNode
			prev.next.next = curr
			self.length += 1

	def sort(self):
		self.head = sortHelper(self.head, self.length - self.length//2, self.length//2)

	def getVal(self, idx):
		curr = self.head
		while (nstep == idx > 0):
			curr = curr.next
		return curr.data

# simple test cases

###################################################
def TestCase1():
	a = SLinkedList()

	print ('\n1')
	a.print()
	print (a.isEmpty())

	print ('\n2')
	a.addHead(4)
	a.addHead(6)
	a.addHead(1)
	a.addEnd(9)
	a.addHead(17)
	a.addEnd(14)
	a.print()
	print (a.isEmpty())
	print (a.length)

	print ('\n3')
	a.insert(3, 12)
	a.insert(2, 23)
	a.insert(0, 19)
	a.addEnd(4)
	a.print()
	print (a.length)

	print ('\n4')
	a.sort()
	a.print()
	print (a.length)

###################################################
if __name__ == "__main__":
	TestCase1()
