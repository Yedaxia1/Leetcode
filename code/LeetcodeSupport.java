import java.util.List;

public class LeetcodeSupport {
}


class TreeNode {
    int val;
    TreeNode left;
    TreeNode right;
    TreeNode() {}
    TreeNode(int val) { this.val = val; }
    TreeNode(int val, TreeNode left, TreeNode right) {
        this.val = val;
        this.left = left;
        this.right = right;
    }
}

class Employee {
    public int id;
    public int importance;
    public List<Integer> subordinates;
};

