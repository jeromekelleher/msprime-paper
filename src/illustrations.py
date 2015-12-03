"""
Python code to generate the illustations in the paper.
"""
from __future__ import print_function
from __future__ import division

import string
import tempfile
import collections
import subprocess

import _msprime
import msprime

def intervals_intersect(a1, b1, a2, b2):
    """
    Returns True if the specified half-closed intervals [a1, b1)
    and [a2, b2) intersect.
    """
    # Intervals do not interset if [a1, b1) is wholly to the left
    # or right of [a2, b2).
    return not (b1 <= a2 or a1 >= b2)

class Tree(object):
    """
    Class representing a simple tree consisting of a set of node->children
    assignments.
    """
    def __init__(self, picture_name, children, parent):
        self.picture_name = picture_name
        self.children = children
        self.parent = parent
        self.num_leaves = len(children) + 1
        self.root = max(children.keys())
        self.x_coords = {}
        self.leaf_x = 0
        self.assign_coordinates(self.root)

    def _get_nodes_to_draw(self, active_nodes):
        if active_nodes is None or len(active_nodes) == 0:
            draw_nodes = set(
                list(self.children.keys()) +
                list(range(1, self.num_leaves + 1))
            )
        else:
            # Do a depth first traversl for all nodes in active_nodes
            draw_nodes = set()
            for u in active_nodes:
                stack = [u]
                while len(stack) > 0:
                    v = stack.pop()
                    draw_nodes.add(v)
                    if v in self.children:
                        for c in self.children[v]:
                            stack.append(c)
        return draw_nodes

    def write_asy_picture(self, out, prefix=None, active_nodes=None,
            x_scale=1, y_scale=1, line_styles=None):
        """
        Writes out the asymptote description of this tree to the
        specified file as a picture.
        """
        ls = collections.defaultdict(lambda: "solid")
        if line_styles is not None:
            ls.update(line_styles)
        draw_nodes = self._get_nodes_to_draw(active_nodes)
        picture_name = self.picture_name
        if prefix is not None:
            picture_name = "{}_{}".format(prefix, picture_name)
        print("picture {};".format(picture_name), file=out)
        print("unitsize({}, 1cm);".format(picture_name), file=out)
        coords = {}
        for u, x in self.x_coords.items():
            y = 0
            if u > self.num_leaves:
                y = u - self.num_leaves
            coords[u] = x * x_scale, y * y_scale
        for u, children in self.children.items():
            for c in children:
                if u in draw_nodes:
                    s = "draw({}, {}--{}, {});".format(
                        picture_name, coords[u], coords[c],
                        ls[u])
                    print(s, file=out)
        for u, z in coords.items():
            if u in draw_nodes:
                s = (
                    "draw({}, '${}$', {}, e=roundbox, "
                    "FillDraw(node_{}));").format(
                        picture_name, u, z, u)
                print(s, file=out)

    def assign_coordinates(self, node):
        """
        Assign x coordinates to all nodes underneath this node.
        """
        if node in self.children:
            for c in self.children[node]:
                self.assign_coordinates(c)
            # We now have x coords for both children
            c1 = self.x_coords[self.children[node][0]]
            c2 = self.x_coords[self.children[node][1]]
            a = min(c1, c2)
            b = max(c1, c2)
            self.x_coords[node] = a + (b - a) / 2
        else:
            self.x_coords[node] = self.leaf_x
            self.leaf_x += 1



class Illustrator(object):
    """
    Superclass of illustrations for coalescent algorithm.
    """
    def __init__(
            self, sample_size, num_loci=10, scaled_recombination_rate=0.2,
            random_seed=1):
        self.tree_sequence = msprime.simulate(
            sample_size, num_loci, scaled_recombination_rate,
            random_seed=random_seed)
        self.sample_size = sample_size
        self.num_loci = num_loci
        self.trees = []
        # first we create the trees.
        C = {}
        pi = {}
        left = 0
        for l, rout, rin in self.tree_sequence.diffs():
            for node, children, _ in rout:
                for c in children:
                    del pi[c]
                del C[node]
            for node, children, _ in rin:
                for c in children:
                    pi[c] = node
                C[node] = children
            picture_name = "tree{}".format(left)
            t = Tree(picture_name, dict(C), dict(pi))
            t.left = left
            t.right = left + l
            self.trees.append(t)
            left += l

    def _generate_pens(self, out):
        """
        Write out the definition of the pens we use to assign
        colours to nodes.
        """
        # Got these from http://tools.medialab.sciences-po.fr/iwanthue/
        # Use the 'intense' preset with k = 9.
        colours = [
                "875373",
                "6AA944",
                "C54C3C",
                "B061D0",
                "4A907F",
                "675927",
                "6F84BD",
                "C8913B",
                "CD4B8A"]
        for j in range(1, self.tree_sequence.get_num_nodes() + 1):
            s = "pen node_{} = rgb(\"{}\");".format(
                    j, colours[j % len(colours)])
            print(s, file=out)

    def write(self, filename, echo=False):
        with tempfile.NamedTemporaryFile() as f:
            self._generate_source(f)
            f.flush()
            f.seek(0)
            if echo:
                print(f.read())
            subprocess.check_call(
                ["asy", "-f", "pdf", "-o", filename, f.name])

    def _generate_source(self, out):
        print("size(0,0);", file=out)
        print("defaultpen(fontsize(8pt));", file=out)
        print("unitsize(1cm);", file=out)
        self._generate_pens(out)

class TreeSequenceIllustrator(Illustrator):
    """
    Illustrates the TreeSequence via an example.
    """
    def __init__(
            self, sample_size, num_loci=10, scaled_recombination_rate=0.2,
            random_seed=1):
        super(TreeSequenceIllustrator, self).__init__(
            sample_size, num_loci, scaled_recombination_rate, random_seed)
        self.fig_height = 5
        self.fig_width = 10
        # Now make the scales.
        self.record_box_height = 0.5 * self.fig_height
        self.tree_box_height = 0.5 * self.fig_height
        times = set()
        for r in self.tree_sequence.records():
            times.add(r[-1])
        self.num_events = len(times)
        max_node = max((t.root) for t in self.trees) - sample_size
        self.fig_x_padding = 0.5
        self.fig_y_padding = 0.1
        self.record_y_scale = self.record_box_height #- (2 * self.fig_y_padding)
        self.record_y_scale /= self.num_events
        self.record_x_scale = self.fig_width - (2 * self.fig_x_padding)
        self.record_x_scale /= num_loci
        self.tree_width = self.fig_width / 5
        self.tree_height = self.tree_box_height - 2 * self.fig_y_padding
        self.tree_x_scale = self.tree_width / sample_size
        self.tree_y_scale = self.tree_height / max_node

    def _get_sparse_tree_label(self, tree):
        max_node = max(t.root for t in self.trees)
        pi = [0 for j in range(max_node + 1)]
        for u, v in tree.parent.items():
            pi[u] = v
        return "".join(str(u) for u in pi[1:])


    def _generate_source(self, out):
        super(TreeSequenceIllustrator, self)._generate_source(out)

        # Draw the genomic scale. TODO this should be generalised
        # and shared with the other figure.
        y_line = 0.25
        tick = 0.05
        x1 = self.fig_x_padding
        x2 = self.fig_width - self.fig_x_padding
        s = "draw({}--{});".format((x1, y_line), (x2, y_line))
        print(s, file=out)
        for l in range(self.num_loci + 1):
            x = x1 + l * self.record_x_scale
            t = tick
            if l % 5 == 0:
                t *= 2
            y1 = y_line - t
            y2 = y_line + t
            s = "draw({}--{});".format((x, y1), (x, y2))
            print(s, file=out)
            if l % 5 == 0:
                s = "label('${}$', {}, N);".format(l, (x, y2))
                print(s, file=out)

        # Print out the trees.
        bottom_line = -self.fig_height + self.fig_y_padding
        for tree in self.trees:
            tree.write_asy_picture(
                out, x_scale=self.tree_x_scale, y_scale=self.tree_y_scale)
            if tree.left != 0:
                x = self.fig_x_padding + self.record_x_scale * tree.left
                y1 = 0
                y2 = -self.fig_height
                s = "draw({}--{}, dashed);".format((x, y1), (x, y2))
                print(s, file=out)
            x = tree.left + (tree.right - tree.left) / 2
            x = self.fig_x_padding + x * self.record_x_scale
            pi = self._get_sparse_tree_label(tree)
            s = "label('$\\\\pi = {}$', {}, S);".format(pi, (x, bottom_line - 0.2))
            print(s, file=out)
            # This _should_ be 2, but we're messing up the scaling somewhere...
            x -= self.tree_width / 2.5
            s = "add(currentpicture, {}.fit(), ({},{}));".format(
                tree.picture_name, x, bottom_line)
            print(s, file=out)

        all_records = list(self.tree_sequence.records())
        d = collections.defaultdict(list)
        for l, r, node, children, time in all_records:
            d[time].append((l, r, node, children))

        # Now write out the records
        records = sorted(list(d.items()), key=lambda x: x[0])
        top_line = -2 * self.fig_y_padding
        x_line = 0.125
        for j, (time, intervals) in enumerate(records):
            time_label = "{:.3f}".format(time)
            picture_name = "currentpicture"
            y = top_line - j * self.record_y_scale
            s = "label({}, '${}$', ({}, {}), W);".format(
                picture_name, time_label, x_line, y)
            print(s, file=out)
            x1 = x_line - tick
            x2 = x_line + tick
            s = "draw({}--{});".format((x1, y), (x2, y))
            print(s, file=out)
            for l, r, node, children in intervals:
                lx = self.fig_x_padding + l * self.record_x_scale
                rx = self.fig_x_padding + r * self.record_x_scale
                s = "draw({}, ({},{})--({},{}));".format(
                    picture_name, lx, y, rx, y)
                print(s, file=out)
                print(
                    "dot({}, ({},{}));".format(picture_name, lx, y),
                    file=out)
                print(
                    "dot({}, ({},{}), filltype=FillDraw(white));".format(
                        picture_name, rx, y),
                    file=out)
                mx = lx + (rx - lx) / 2
                s = "label({}, '${}\\\\rightarrow {}$', ({}, {}), N);".format(
                        picture_name, children, node, mx, y)
                print(s, file=out)


        # Draw the time scale line
        y1 = top_line
        y2 = top_line - (self.num_events - 1) * self.record_y_scale
        s = "draw({}--{});".format((x_line, y1), (x_line, y2))
        print(s, file=out)


class AlgorithmIllustrator(Illustrator):
    """
    Class producing illustrations of Hudson's algorithm.
    """
    def __init__(
            self, sample_size, num_loci=10, scaled_recombination_rate=0.2,
            random_seed=1):
        super(AlgorithmIllustrator, self).__init__(
            sample_size, num_loci, scaled_recombination_rate, random_seed)
        self.simulator = _msprime.Simulator(
            sample_size=sample_size, num_loci=num_loci,
            scaled_recombination_rate=scaled_recombination_rate,
            random_seed=random_seed)
        self.box_height = 3.8
        self.box_width = 6
        self.num_cols = 2
        self._run_simulation()
        # Now make the scales.
        self.state_box_height = 0.5 * self.box_height
        self.tree_box_height = 0.6 * self.box_height
        # We need to know the maximum number of ancestors and tree nodes
        # so we can make the appropriate scales.
        max_ancestors = max(len(a) for a in self.states)
        max_node = max((t.root) for t in self.trees) - sample_size
        self.box_x_padding = 0.5
        self.box_y_padding = 0.25
        self.ancestor_y_scale = self.state_box_height - (2 * self.box_y_padding)
        self.ancestor_y_scale /= max_ancestors
        self.ancestor_x_scale = self.box_width - (2 * self.box_x_padding)
        self.ancestor_x_scale /= num_loci
        self.tree_width = self.box_width / 6
        self.tree_height = self.tree_box_height - 2 * self.box_y_padding
        self.tree_x_scale = self.tree_width / sample_size
        self.tree_y_scale = self.tree_height / max_node

    def _run_simulation(self):
        self.states = []
        self.breakpoints = []
        self.event_types = []
        self.times = []
        self.ancestor_ids = {}
        self.next_ancestor_id = 1
        sim = self.simulator
        sim.run(-1)
        not_done = True
        event_type = ""
        while not_done:
            ancestors = sim.get_ancestors()
            self.states.append(ancestors)
            self.breakpoints.append(
                [0] + self.simulator.get_breakpoints() + [self.num_loci])
            self.times.append(sim.get_time())
            self.event_types.append(event_type)
            before = sim.get_num_recombination_events()
            not_done = not sim.run_event()
            event_type = "CA"
            if before != sim.get_num_recombination_events():
                event_type = "RE"
        self.states.append(sim.get_ancestors())
        self.breakpoints.append(
            [0] + self.simulator.get_breakpoints() + [self.num_loci])
        self.times.append(sim.get_time())
        self.event_types.append(event_type)

    def _get_intervals(self, ancestors, breakpoints):
        """
        Returns the set of left, right intervals in the set along with the
        nodes present in this interval.
        """
        intervals = collections.OrderedDict()
        for j in range(len(breakpoints) - 1):
            intervals[breakpoints[j], breakpoints[j + 1]] = []
        intervals[breakpoints[-1], self.num_loci] = []
        # Now we collect the nodes that are represented in each interval.
        for (a1, b1), l in intervals.items():
            for ancestor in ancestors:
                for a2, b2, node in ancestor:
                    if intervals_intersect(a1, b1, a2, b2):
                        l.append(node)
        return intervals

    def _get_num_links(self, ancestors):
        """
        Returns the number of recombination links.
        """
        num_links = 0
        for ancestor in ancestors:
            l = ancestor[0][0]
            r = ancestor[-1][1]
            num_links += r - l - 1
        return num_links

    def _get_ancestor_id(self, ancestor):
        a = tuple(ancestor)
        if a not in self.ancestor_ids:
            self.ancestor_ids[a] = self.next_ancestor_id
            self.next_ancestor_id += 1
        return self.ancestor_ids[a]

    def _generate_state_picture(
            self, picture_name, ancestors, breakpoints, event_type,
            time, label, out):
        """
        Generates the picture representing the state for the
        specified set of ancestors.
        """
        print("picture {};".format(picture_name), file=out)
        print("unitsize({}, 1cm);".format(picture_name), file=out)
        print("size({}, {}, {});".format(
            picture_name, self.box_width, self.box_height), file=out)
        x = self.box_x_padding
        y = -self.box_y_padding
        num_links = self._get_num_links(ancestors)
        s = "label({}, '$({})$ $t={:.3f} \qquad L={}$', ({}, {}), E);".format(
                picture_name, label, time, num_links, x, y)
        x = self.box_width - self.box_x_padding
        print(s, file=out)
        s = "label({}, '{}', ({}, {}), W);".format(
                picture_name, event_type, x, y)
        print(s, file=out)
        top_line = -2.5 * self.box_y_padding
        for j, ancestor in enumerate(reversed(ancestors)):
            y = top_line - j * self.ancestor_y_scale
            # First draw the dotted line to indicate the full genomic range.
            lx = self.box_x_padding
            rx = self.box_x_padding +  self.num_loci * self.ancestor_x_scale
            s = "draw({}, ({},{})--({},{}), dotted);".format(
                picture_name, lx, y, rx, y)
            print(s, file=out)
            s = "label({}, '$l_{{{}}}$', ({}, {}), W);".format(
                picture_name, self._get_ancestor_id(ancestor), lx, y)
            print(s, file=out)
            for l, r, node in ancestor:
                lx = self.box_x_padding + l * self.ancestor_x_scale
                rx = self.box_x_padding + r * self.ancestor_x_scale
                s = "draw({}, ({},{})--({},{}), node_{} + 2);".format(
                    picture_name, lx, y, rx, y, node)
                print(s, file=out)
                print(
                    "dot({}, ({},{}));".format(picture_name, lx, y),
                    file=out)
                print(
                    "dot({}, ({},{}), filltype=FillDraw(white));".format(
                        picture_name, rx, y),
                    file=out)
                node_label = "{}".format(node)
                s = (
                    "draw({}, '${}$', ({}, {}), e=roundbox, "
                    "FillDraw(node_{}));".format(
                        picture_name, node_label, lx + (rx - lx) / 2, y, node)
                )
                print(s, file=out)

        # Print out the partial trees.
        bottom_line = -self.box_height + self.box_y_padding
        for interval, active_nodes in self._get_intervals(
                ancestors, breakpoints).items():
            for tree in self.trees:
                if tree.left == interval[0]:
                    tree.write_asy_picture(
                        out, picture_name, active_nodes, self.tree_x_scale,
                        self.tree_y_scale)
                    name = "{}_{}".format(picture_name, tree.picture_name)
                    x = interval[0] + (interval[1] - interval[0]) / 2
                    x = self.box_x_padding + x * self.ancestor_x_scale
                    s = "add({}, {}.fit(), ({},{}), N);".format(
                        picture_name, name, x, bottom_line)
                    print(s, file=out)

    def _generate_source(self, out):
        super(AlgorithmIllustrator, self)._generate_source(out)
        num_states = len(self.states)
        # Draw the grid around the states.
        num_cols = self.num_cols
        num_rows = num_states // self.num_cols
        for j in range(num_rows + 1):
            y = -j * self.box_height
            s = "draw({}--{});".format((0, y), (num_cols * self.box_width, y))
            print(s, file=out)
        for j in range(num_cols + 1):
            x = j * self.box_width
            s = "draw({}--{});".format((x, 0), (x, -num_rows * self.box_height))
            print(s, file=out)
        # Draw the genomic scales.
        y_line = 0.25
        tick = 0.05
        for j in range(num_cols):
            x1 = self.box_x_padding + self.box_width * j
            x2 = self.box_width * (j + 1) - self.box_x_padding
            s = "draw({}--{});".format((x1, y_line), (x2, y_line))
            print(s, file=out)
            for l in range(self.num_loci + 1):
                x = x1 + l * self.ancestor_x_scale
                t = tick
                if l % 5 == 0:
                    t *= 2
                y1 = y_line - t
                y2 = y_line + t
                s = "draw({}--{});".format((x, y1), (x, y2))
                print(s, file=out)
                if l % 5 == 0:
                    s = "label('${}$', {}, N);".format(l, (x, y2))
                    print(s, file=out)
        letters = string.ascii_lowercase
        for j in range(num_states):
            ancestors = self.states[j]
            breakpoints = self.breakpoints[j]
            event_type = self.event_types[j]
            time = self.times[j]
            picture_name = "state{}".format(j)
            self._generate_state_picture(
                picture_name, ancestors, breakpoints, event_type,
                time, letters[j], out)
            # Now we must place the picture.
            y = -(j // num_cols) * self.box_height
            x = (j % num_cols) * self.box_width
            s = "add(currentpicture, {}.fit(), ({},{}));".format(
                picture_name, x, y)
            print(s, file=out)


def main():

    illustrator = TreeSequenceIllustrator(4, random_seed=21)
    illustrator.write("figures/tree-sequence-illustration.pdf", False)

    illustrator = AlgorithmIllustrator(4, random_seed=21)
    illustrator.write("figures/hudsons-algorithm-illustration.pdf", False)

if __name__ == "__main__":
    main()
