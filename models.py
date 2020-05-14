

class linkedNode:
    """
    linked node object that keeps track of the value of the node itself, the value of the source node,
    and the direction of the source node
    """
    def __init__(self, value, source, direction):
        self.value = value
        self.source = source
        self.direction = direction
