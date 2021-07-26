import warnings

class DEDlattice:
    def __init__(self, dim):
        self.dim = dim
        facfile = f'data/m{self.dim}.fac'
        hpfile = f'data/m{self.dim}.out'
        self.ps = KunzPoset.ReadFacesFromNormaliz(face_lattice_file_path=facfile, hplane_file_path=hpfile)
        self.pd = self.classify()

    def get_ps(self):
        return self.ps

    def get_pd(self):
        return self.pd

    def check_refinement(self, p1, p2):
        """Given two posets p1 and p2, determine if p1 contains all relations in P2. If so, P1 is a refinement of P2.

        Input can be poset objects, or int values of position in ps.

        :param p1: poset p1
        :type p1: int or poset object
        :param p2: poset p2
        :type p2: int or poset object
        :return: if p1 is a refinement of p2
        :rtype: bool
        """
        if isinstance(p1, sage.rings.integer.Integer) and isinstance(p2, sage.rings.integer.Integer):
            r1 = self.ps[p1].hyperplane_desc
            r2 = self.ps[p2].hyperplane_desc
        else:
            r1 = p1.hyperplane_desc
            r2 = p2.hyperplane_desc
        for r in r2:
            if r not in r1:
                return False
        return True

    def classify(self):
        """Classify all posets of the Kunz polyhedron by their dimensions.

        :return: A list of the posets ordered by dimension (i-th list corresponds to dimension i).
        :rtype: list of lists
        """
        result = []
        for i in range(self.ps[0].m):
            result.append([])

        for i in range(len(self.ps)):
            result[self.ps[i].Dimension()].append(i)

        return result

    def child_posets(self, p):
        """Find all posets below a given poset p.

        Input can be poset object or int value of its position in ps.

        :param p: poset p
        :type p: int or poset object
        :return: list of integers, representing positions of corresponding posets
        :rtype: list
        """
        if isinstance(p, sage.rings.integer.Integer):
            p = self.ps[p]
        result = []
        d = p.Dimension()
        for i in self.pd[d - 1]:
            if self.check_refinement(self.ps[i], p):
                result.append(i)
        return result

    def child_ed(self, p):
        """For a given poset p, find all posets below p, including their embedding dimension
        and number of defining facet equations.

        Input can be poset object or int value of its position in ps.

        :param p: poset p
        :type p: int or poset object
        :return: list of poset indices, embedding dimensions and #(def. facet equations)
        :rtype: list of 3-tuples (i, e, n)
        """
        if isinstance(p, sage.rings.integer.Integer):
            p = self.ps[p]
        result = []
        d = p.Dimension()
        for i in self.pd[d - 1]:
            if self.check_refinement(self.ps[i], p):
                result.append((i, len(self.ps[i].atoms), len(self.ps[i].hyperplane_desc)))
        return result

    def child_relations(self, p):
        """For a given poset p, find the number of facets added for all its child posets.

        Input can be poset object or int value of its position in ps.

        :param p: poset p
        :type p: int or poset object
        :return: list of poset indices and #added facet equations
        :rtype: list of 2-tuples (i,n)
        """
        if isinstance(p, sage.rings.integer.Integer):
            p = self.ps[p]
        result = []
        d = p.Dimension()
        n = len(p.hyperplane_desc)      # Is this supposed to do something?
        for i in self.pd[d - 1]:
            if self.check_refinement(self.ps[i], p):
                result.append((i, len(self.ps[i].hyperplane_desc)))
        return result

    def ed_poset(self):
        """Compute the cover relations for the DED lattice of the Kunz polyhedron.

        :return: Cover relations of the DED lattice, stored as a list of ((d1,e1),(d2,e2)),
                    where d1 is the face dim and e1 the embedding dim of the child poset
                    and d2 is the face dim and e2 the embedding dim of the parent poset.
        :rtype: list of 2-tuples
        """
        result = []
        for i in range(len(self.ps)):
            p = self.ps[i]
            d = p.Dimension()
            ed = len(p.atoms)
            ef = len(p.hyperplane_desc)     # Is this supposed to do something?
            lchild = self.child_ed(p)
            for (c, cd, cf) in lchild:
                if ((d - 1, cd), (d, ed)) not in result:
                    result.append(((d - 1, cd), (d, ed)))
        return result

    def ef_poset(self):
        """Compute the cover relations for the DEF lattice of the Kunz polyhedron.

        :return: Cover relations of the DEF lattice, stored as a list of ((d1,e1,n1),(d2,e2,n2)),
                    where d1 is the face dim, e1 the embedding dim, n1 the #(def. facet equations) of the child poset
                    and d2 is the face dim, e2 the embedding dim, n2 the #(def. facet equations) of the parent poset.
        :rtype: list of 2-tuples
        """
        result = []
        for i in range(len(self.ps)):
            p = self.ps[i]
            d = p.Dimension()
            ed = len(p.atoms)
            ef = len(p.hyperplane_desc)
            lchild = self.child_ed(p)
            for (c, cd, cf) in lchild:
                if ((d - 1, cd, cf), (d, ed, ef)) not in result:
                    result.append(((d - 1, cd, cf), (d, ed, ef)))
        return result

    def relation_poset(self):
        """Compute the cover relations for the simplified DEF lattice of the Kunz polyhedron,
        without embedding dimension.

        :return: Cover relations of the DEF lattice, stored as a list of ((d1,n1),(d2,n2)),
                    where d1 is the face dim, n1 the #(def. facet equations) of the child poset
                    and d2 is the face dim, n2 the #(def. facet equations) of the parent poset.
        :rtype: list of 2-tuples
        """
        result = []
        for i in range(len(self.ps)):
            p = self.ps[i]
            d = p.Dimension()
            ef = len(p.hyperplane_desc)
            lchild = self.child_relations(p)
            for (c, cf) in lchild:
                if ((d - 1, cf), (d, ef)) not in result:
                    result.append(((d - 1, cf), (d, ef)))
        return result

    def find_poset_nd(self, d, n):
        """Find posets with n defining facet equations and face dimension d.

        :param d: face dimension
        :type d: int
        :param n: #(def. facet equations)
        :type n: int
        :return: indices of posets
        :rtype: list
        """
        result = []
        pd_d = self.pd[d]
        for i in pd_d:
            if len(self.ps[i].hyperplane_desc) == n:
                result.append(i)
        return result

    def find_poset_ed(self, d, e):
        """Find posets with face dimension d and embedding dimension e.

        :param d: face dimension
        :type d: int
        :param e: embedding dimension
        :type e: int
        :return: indices of posets with embedding dim e and face dim d
        :rtype: list
        """
        result = []
        pd_d = self.pd[d]
        for i in pd_d:
            if len(self.ps[i].atoms) == e:
                result.append(i)
        return result

    def add_one_facet(self, p):
        """Find a poset below p by adding a single facet equation, if possible.

        :param p: poset p
        :type p: int or poset object
        :return: the index of the found poset, or False if not possible
        :rtype: int or bool
        """
        if isinstance(p, sage.rings.integer.Integer):
            p = self.ps[p]
            cd = self.child_ed(p)
            for i in range(len(cd)):
                if cd[i][2] - (len(p.hyperplane_desc)) == 1:
                    return cd[i][0]
            return False

    # this function is *supposed* to find a flag that contains the maximal jump of embedding dimension
    # (in consturction or unnecessary for now)
    def find_maximal_flag(self):
        warnings.warn("WARNING: This method is under construction.")
        result = [(0, self.ps[0].m - 1, self.ps[0].m - 1)]
        p = self.ps[0]
        while len(result) != self.ps[0].m - 2:
            cd = self.child_ed(p)
            # maxd = -1
            # c = 0
            # for i in range(len(cd)):
            #    if cd[i][1] - (Ps[0].m-len(result)) > maxd:
            #        c = cd[i][0]
            #        maxd = cd[i][1] - (Ps[0].m-len(result))
            # result.append((c, Ps[0].m-len(result), cd[i][1]))
            # P = Ps[c]
            for i in range(len(cd)):
                if cd[i][2] - (len(result) - 1) == 1:
                    p = self.ps[cd[i][0]]
                    result.append((cd[i][0], self.ps[0].m - len(result) - 1, cd[i][1]))
                    break
                if i == len(cd) - 1:
                    return result

        cd = self.child_ed(p)
        for i in range(len(cd)):
            if cd[i][1] == 1:
                result.append((cd[i][0], 1, 1))
                break
        return result