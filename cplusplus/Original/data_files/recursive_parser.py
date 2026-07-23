import re

class Parser:
    def __init__(self, text):
        # Tokenize parentheses and words
        self.tokens = re.findall(r'\(|\)|[A-Za-z0-9_]+', text)
        self.pos = 0

    def peek(self):
        return self.tokens[self.pos]

    def consume(self):
        tok = self.tokens[self.pos]
        self.pos += 1
        return tok

    def parse(self):
        tok = self.consume()

        # Variable
        if tok != "(":
            return tok

        gate = self.consume()

        if gate == "WIRE":
            expr = self.parse()
            self.consume()  # )
            return expr

        elif gate == "NOT1":
            expr = self.parse()
            self.consume()
            return f"¬({expr})"

        elif gate in ("AND1", "OR1", "XOR1"):
            left = self.parse()
            right = self.parse()
            self.consume()

            op = {
                "AND1": "*",
                "OR1": "+",
                "XOR1": "⊕"
            }[gate]

            return f"({left} {op} {right})"

        else:
            raise ValueError(f"Unknown gate: {gate}")


if __name__ == "__main__":
    expression = "(AND1 (XOR1 (OR1 (WIRE (AND1 A B)) (OR1 (XOR1 A C) (XOR1 A C))) (OR1 (NOT1 (OR1 B A)) (WIRE (WIRE D)))) (XOR1 (OR1 (WIRE (AND1 A B)) (OR1 (XOR1 A C) (XOR1 A C))) (OR1 (NOT1 (OR1 B A)) (WIRE (WIRE D)))))"

    parser = Parser(expression)
    formula = parser.parse()

    print(formula)