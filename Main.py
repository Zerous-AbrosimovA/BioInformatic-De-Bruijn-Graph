from collections import defaultdict
from typing import Any
from Bio import SeqIO

class BuildGraph:
    def __init__(self, k: int, outputFASTA: str, outputGFA: str) -> None:
        self.k = k
        self.outputFASTA = outputFASTA
        self.outputGFA = outputGFA
        self.graph = defaultdict(lambda: defaultdict(lambda: ""))
        self.depth = defaultdict(lambda: defaultdict(lambda: 0))
        self.inDegree = defaultdict(lambda: 0)
        self.outDegree = defaultdict(lambda: 0)
        self.visited = set()

    def add(self, reading: str) -> None:
        for i in range(len(reading) - self.k):
            u: str = reading[i : i + self.k + 1][:-1]
            v: str = reading[i : i + self.k + 1][1:]
            if v not in self.graph[u]:
                self.inDegree[v] += 1
                self.outDegree[u] += 1
                _ = self.graph[v]
            self.graph[u][v] = str(reading[i : i + self.k + 1])
            self.depth[u][v] += 1

    def createCompressedGraph(self) -> None:
        self._compressGraph()

    def createClearedGraph(self) -> None:
        self._clearGraph()

    def printGraph(self) -> None:
        self._printFASTA()
        self._printGFA()

    def _compressGraph(self) -> None:
        vertexes = self.graph.keys()
        compressedGraph: list = []
        startingPoints: list = [i for i in vertexes if self.inDegree[i] != 1 or self.outDegree[i] != 1]
        for u in startingPoints:
            for v in self.graph[u]:
                if not self.visited.__contains__((u, v)):
                    compressedGraph.append(self._iterativeDfs(u, v))
        for u in self.graph.keys():
            for v in self.graph[u]:
                if not self.visited.__contains__((u, v)):
                    compressedGraph.append(self._iterativeDfs(u, v))
        self._rebuildGraph(compressedGraph)


    def _iterativeDfs(self, u: str, v: str) -> dict[str, str | float | int | Any]:
        start = u
        currentPath: str = self.graph[u][v]
        depth: int = self.depth[u][v]
        totalEdges: int = 1
        self.visited.add((u, v))
        while (self.inDegree[v] == 1 and self.outDegree[v] == 1 and
               not self.visited.__contains__((v, list(self.graph[v].keys())[0]))):
            neighbor: str = list(self.graph[v].keys())[0]
            self.visited.add((v, neighbor))
            currentPath += self.graph[v][neighbor][-1]
            depth += self.depth[v][neighbor]
            totalEdges += 1
            u = v
            v = neighbor
        return {"start" : start, "end" : v, "nucleotides" : currentPath,
                "coverage" : depth / totalEdges, "length" : totalEdges}

    def _clearGraph(self) -> None:
        innerGraph = []
        for u in self.graph:
            for v in self.graph[u]:
                innerGraph.append({"start": u, "end": v,
                    "nucleotides": self.graph[u][v], "coverage": self.depth[u][v], "length": 0})
        totalCoverage: int = 0
        for i in innerGraph:
            totalCoverage += i["coverage"]
        avgCoverage: float = totalCoverage / len(innerGraph)
        clearedGraph: list = []
        for i in innerGraph:
            if i["coverage"] >= avgCoverage * 0.3 and (self.inDegree[i["start"]] != 0 or self.outDegree[i["end"]] != 0):
                clearedGraph.append(i)
        self._rebuildGraph(clearedGraph)

    def _rebuildGraph(self, currentGraph) -> None:
        self.graph = defaultdict(lambda: defaultdict(lambda: ""))
        self.depth = defaultdict(lambda: defaultdict(lambda: 0))
        self.inDegree = defaultdict(lambda: 0)
        self.outDegree = defaultdict(lambda: 0)
        self.visited = set()
        for i in currentGraph:
            u = i["start"]
            v = i["end"]
            self.graph[u][v] = i["nucleotides"]
            self.depth[u][v] = i["coverage"] * i["length"]
            self.inDegree[v] += 1
            self.outDegree[u] += 1

    def _printFASTA(self) -> None:
        index: int = 1
        with open(self.outputFASTA, "w") as f:
            for i in self.graph:
                for j in self.graph[i]:
                    f.write(f">index {index} NC_000913.3 Escherichia coli str. K-12 substr. MG1655\n")
                    f.write(self.graph[i][j] + "\n")
                    index += 1

    def _printGFA(self) -> None:
        with open(self.outputGFA, "w") as f:
            f.write("H\tVN:Z:1.1\n")
            starting_at = defaultdict(list)
            index: int = 0
            for i in self.graph:
                for j in self.graph[i]:
                    f.write(f"S\t{index}\t{self.graph[i][j]}\n")
                    starting_at[i].append(index)
                    index += 1
            index = 0
            for i in self.graph:
                for j in self.graph[i]:
                    if j in starting_at:
                        for neighbor in starting_at[j]:
                            f.write(f"L\t{index}\t+\t{neighbor}\t+\t{self.k - 1}M\n")
                            index += 1


graphClass = BuildGraph(51, "output.fasta", "output.gfa")
for record in SeqIO.parse("ecoli_reads.fastq.txt", "fastq"):
    graphClass.add(record.seq)

graphClass.createCompressedGraph()
graphClass.createClearedGraph()
graphClass.createCompressedGraph()
graphClass.printGraph()
