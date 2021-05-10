import javax.persistence.criteria.CriteriaBuilder;
import java.util.*;

public class Leetcode {
    // # 4
    public double findMedianSortedArrays(int[] nums1, int[] nums2) {
        if(nums1.length==0)
            return nums2.length%2==0? (((double)(nums2[nums2.length/2]+nums2[nums2.length/2-1]))/(double)2):(double)(nums2[nums2.length/2]);
        if(nums2.length==0)
            return nums1.length%2==0? (((double)(nums1[nums1.length/2]+nums1[nums1.length/2-1]))/(double)2):(double)(nums1[nums1.length/2]);

        int i=0, j=0;
        boolean flag = (nums1.length + nums2.length)%2==0? true:false;

        int count = (nums1.length + nums2.length + 1) / 2 - 1;

        if(nums1.length>nums2.length){
            int[] temp = nums1.clone();
            nums1 = nums2.clone();
            nums2 = temp.clone();
        }

        while (count!=0){
            if(i==nums1.length){
                 j += count;
                 return flag? (double)(nums2[j] + nums2[j+1])/(double)(2):(double)nums2[j];
            }

            if(nums1[i] <= nums2[j]){
                i ++;
                count --;
            }else {
                j ++;
                count --;
            }


        }
        if(i==nums1.length){
            j += count;
            return flag? (double)(nums2[j] + nums2[j+1])/(double)(2):(double)nums2[j];
        }
        if(nums1[i] < nums2[j])
            return flag? (double)(nums1[i] + Math.min((i+1 < nums1.length? nums1[i+1]:Integer.MAX_VALUE), nums2[j])) / 2.0 : nums1[i];
        else {
            return flag? (double)(nums2[j] + Math.min(nums1[i], (j+1 < nums2.length? nums2[j+1]:Integer.MAX_VALUE))) / 2.0 : nums2[j];
        }


    }

    // # 5
    public String longestPalindrome(String s) {
        if(s=="" || s.length()==1)
            return s;
        int max_len = 1;
        int ans_index = 0;
        boolean[][] dp = new boolean[s.length()][s.length()];

        for (int i=0;i<s.length();i++)
            dp[i][i] = true;

        for(int l=2;l<=s.length();l++){
            for (int i=0;i<s.length()-l+1;i++){
                dp[i][i+l-1] = (s.charAt(i)==s.charAt(i+l-1)? true:false) && (l>2? dp[i+1][i+l-2]:true);
                if(dp[i][i+l-1]){
                    if(l > max_len){
                        max_len = l;
                        ans_index = i;
                    }
                }
            }
        }

        return s.substring(ans_index, ans_index + max_len);
    }

    // # 6
    public String convert(String s, int numRows) {
        if(s.length()==1)
            return s;
        ArrayList[] lists = new ArrayList[numRows];
        lists[0].add(s.charAt(0));
        for (int i=1;i<s.length();i++){
            lists[(i % (numRows + 1))].add(s.charAt(i));
        }
        StringBuffer stringBuffer = new StringBuffer();
        for (List<Character> list : lists)
            stringBuffer.append(list.toString());
        return stringBuffer.toString();
    }

    // # 7
    public int reverse(int x) {
        if(x <= Integer.MIN_VALUE || x >= Integer.MAX_VALUE || x==0)
            return 0;
        int flag = x < 0? -1:1;
        int ans = 0;
        x = Math.abs(x);
        while (x != 0){
            if(ans > Integer.MAX_VALUE / 10)
                return 0;
            ans = ans * 10 + x%10;
            //System.out.println(ans);
            x /= 10;
        }
        return flag * ans;
    }

    // # 11
    public int maxArea(int[] height) {
        int left = 0, right = height.length-1;
        int ans = Integer.MIN_VALUE;
        while (left < right){
            if(Math.min(height[left], height[right]) * (right-left) > ans)
                ans = Math.min(height[left], height[right]) * (right-left);
            if (height[left] < height[right]) {
                left++;
            } else {
                right--;
            }
        }
        return ans;
    }

    // # 27
    public int removeElement(int[] nums, int val) {
        int ans = 0;
        int left = 0;
        int right = nums.length - 1;
        while (left<right){
            if(nums[left] != val){
                left++;
                continue;
            }else {
                while (nums[right] == val && left<right){
                    right--;
                }
                if(left<right){
                    int temp = nums[left];
                    nums[left] = nums[right];
                    nums[right] = temp;
                }
            }
        }
        return left;
    }

    // # 28
    public int strStr(String haystack, String needle) {
        if(needle.equals(""))
            return 0;
        for(int i=0; i < haystack.length()-needle.length(); i++){
            if(needle.compareTo(haystack.substring(i, i+needle.length()))==0)
                return i;
        }
        return -1;
    }

    // # 80
    public int removeDuplicates(int[] nums) {
        int n = nums.length;
        if(n<=2)
            return n;

        int slow = 2, fast = 2;
        while (fast < n){
            if(nums[fast]!=nums[slow-2]){
                nums[slow] = nums[fast];
                slow++;
            }
            fast++;
        }
        return slow;


    }

    // # 81
    public boolean search(int[] nums, int target) {
        if(nums[0]==target)
            return true;

        if(nums.length==1)
            return nums[0]==target;

        if(nums[0] < target){
            int flag = 1;
            while (true){
                if(nums[flag]==target)
                    return true;
                else if(nums[flag] > target)
                    return false;
                flag++;
                if(flag>=nums.length)
                    return false;
                if(nums[flag] < nums[flag-1])
                    return false;
            }
        }else {
            int flag = nums.length-1;
            while (true){
                if(nums[flag]==target)
                    return true;
                else if(nums[flag] < target)
                    return false;
                flag--;
                if(flag<0)
                    return false;
                if(nums[flag] > nums[flag+1])
                    return false;
            }
        }
    }

    // # 91
    public int numDecodings(String s) {
        if(s.charAt(0)=='0')
            return 0;
        if(s.length()==0)
            return 0;
        if (s.length()==1)
            return 1;
        int[] dp = new int[s.length()];
        dp[0] = 1;
        if(s.charAt(1)=='0'){
            if (s.charAt(0)=='1' || s.charAt(0)=='2')
                dp[1] = 1;
            else return 0;
        }else {
            if ((1 <= Integer.parseInt(s.substring(0, 2)) && Integer.parseInt(s.substring(0, 2)) <= 26))
                dp[1] = 2;
            else dp[1] = 1;
        }
        for (int i=2;i<s.length();i++){
            if(s.charAt(i)=='0') {
                if(1 <= Integer.parseInt(s.substring(i-1, i+1)) && Integer.parseInt(s.substring(i-1, i+1)) <= 26)
                    dp[i] = dp[i-2];
                else
                    return 0;
            }
            else {
                if(s.charAt(i-1)=='0')
                    dp[i] = dp[i-1];
                else {
                    dp[i] = dp[i-1];
                    if(Integer.parseInt(s.substring(i-1, i+1)) <= 26)
                        dp[i] += dp[i-2];
                }
            }

        }
        return dp[s.length()-1];
    }

    // # 137
    public int singleNumber(int[] nums) {
        if(nums.length==1)
            return nums[0];
        Arrays.sort(nums);
        int i = 0;
        while (true){
            if(i < nums.length-3 && nums[i]==nums[i+1] && nums[i]==nums[i+2])
                i = i + 3;
            else return nums[i];
        }
    }

    // # 153
    public int findMin(int[] nums) {
        if(nums.length==1)
            return nums[0];
        int ans = nums[0];

        for(int i=1;i<nums.length;i++){
            ans = Math.min(nums[i], ans);

            if(nums[i-1] > nums[i])
                break;
        }
        return ans;

    }

    // # 154
    public int findMin2(int[] nums) {
        if(nums.length==1)
            return nums[0];
        int ans = nums[0];

        for(int i=1;i<nums.length;i++){
            ans = Math.min(nums[i], ans);

            if(nums[i-1] > nums[i])
                break;
        }
        return ans;
    }

    // # 179
    public String largestNumber(int[] nums) {
        if(nums.length==1)
            return String.valueOf(nums[0]);
        int sum = 0;
        for(int i:nums)
            sum+=i;
        if(sum==0)
            return "0";
        Integer[] integers = new Integer[nums.length];
        for(int i=0;i<nums.length;i++)
            integers[i] = nums[i];

        Arrays.sort( integers, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                String s1 = String.valueOf(o1);
                String s2 = String.valueOf(o2);
                if ((int) s1.charAt(0) > (int) s2.charAt(0))
                    return -1;
                else if ((int) s1.charAt(0) < (int) s2.charAt(0)) {
                    return 1;
                } else {
                    StringBuffer ss1 = new StringBuffer();
                    StringBuffer ss2 = new StringBuffer();
                    ss1.append(s1);
                    ss1.append(s2);
                    ss2.append(s2);
                    ss2.append(s1);
                    for(int k=0;k<s1.length()+s2.length();k++){
                        if((int)ss1.charAt(k) < (int)ss2.charAt(k))
                            return 1;
                        else if((int)ss1.charAt(k) > (int)ss2.charAt(k))
                            return -1;
                        else continue;
                    }
                    return -1;
                }
            }
        });

        StringBuffer ans = new StringBuffer();
        for(int i: integers)
            ans.append(i);

        return ans.toString();

    }

    // # 213
    public int rob(int[] nums) {
        if(nums.length==1)
            return nums[0];
        if(nums.length==2)
            return Math.max(nums[0], nums[1]);
        int ans = 0;
        int[][] dp = new int[2][nums.length+1];
        dp[0][1] = 0;
        dp[1][1] = nums[0];
        for(int i=1;i<nums.length-1;i++){
            dp[0][i+1] = Math.max(dp[1][i], dp[0][i]);
            dp[1][i+1] = dp[0][i] + nums[i];
        }
        ans = Math.max(dp[0][nums.length-1], dp[1][nums.length-1]);

        dp[0][2] = 0;
        dp[1][2] = nums[1];
        for(int i=2;i<nums.length;i++){
            dp[0][i+1] = Math.max(dp[1][i], dp[0][i]);
            dp[1][i+1] = dp[0][i] + nums[i];
        }
        ans = Math.max(Math.max(dp[0][nums.length], dp[1][nums.length]), ans);
        return ans;
    }

    // # 363
    public int maxSumSubmatrix(int[][] matrix, int k) {
        int max_ls_k = Integer.MIN_VALUE;
        int[][] line_sum = new int[matrix.length+1][matrix[0].length];
        for (int i=0;i<matrix[0].length;i++)
            line_sum[0][i] = 0;

        for(int i=1;i<matrix.length+1;i++)
            for (int j=0;j<matrix[0].length;j++)
                line_sum[i][j] = matrix[i-1][j] + line_sum[i-1][j];
        // printMatr(line_sum);
        for (int i=1;i<matrix.length+1;i++){
            for (int j=i;j<matrix.length+1;j++){
                int[] sub_sum = new int[matrix[0].length];
                for(int t=0;t<matrix[0].length;t++)
                    sub_sum[t] = line_sum[j][t] - line_sum[i-1][t];
                // System.out.println(Arrays.toString(sub_sum));
                int[] sum_index = new int[sub_sum.length+1];
                sum_index[0] = 0;
                for (int t=1;t<sub_sum.length+1;t++)
                    sum_index[t] = sum_index[t-1] + sub_sum[t-1];
                // System.out.println(Arrays.toString(sum_index));
                for (int m=1;m<sub_sum.length+1;m++)
                    for (int n=m;n<sub_sum.length+1;n++){
                        if((sum_index[n] - sum_index[m-1]) <=k && (sum_index[n] - sum_index[m-1]) > max_ls_k)
                            max_ls_k = sum_index[n] - sum_index[m-1];
                    }
            }
        }
        return max_ls_k;
    }

    // # 373
    public int combinationSum4(int[] nums, int target) {
        // Arrays.sort(nums);
        int[] dp = new int[target+1];
        for (int i=0;i<target+1;i++){
            dp[i] = 0;
        }
        dp[0] = 1;
        for (int i=0;i<target+1 ;i++){
            for (int n:nums){
                if(n<=i){
                    dp[i] += dp[i-n];
                }
            }
        }
        return dp[target];
    }

    // # 403
    public boolean canCross(int[] stones) {
        if(stones[1]!=1)
            return false;
        return jump(0, 1, stones);
    }

    public boolean jump(int now, int step, int[] stones){
        if(now + step == stones[stones.length-1])
            return true;
        if(now + step < stones[stones.length-1])
            if(contain_num(now + step, stones))
                return jump(now+step, step-1, stones)
                        || jump(now+step, step, stones)
                        || jump(now+step, step+1, stones);
        return false;
    }

    public boolean contain_num(int n, int[] stones){
        int left = 0;
        int right = stones.length-1;
        while (left < right){

            int mid = (left + right) / 2;
            System.out.println(left + "--" + mid + "--" + right);
            if(stones[mid]==n)
                return true;
            if(stones[mid] < n)
                left = mid + 1;
            else right = mid - 1;
        }
        return stones[right]==n? true:false;
    }

    // # 554
    public int leastBricks(List<List<Integer>> wall) {
        int hangshu = wall.size();
        int lieshu = wall.get(0).stream().reduce(Integer::sum).orElse(0);
        HashMap<Integer, Integer> hashMap =  new LinkedHashMap<>();
        for (List<Integer> list:wall){
            int nowindex = 0;
            for (Integer integer:list){
                nowindex += integer;
                if(nowindex==lieshu)
                    break;
                if(!hashMap.containsKey(nowindex))
                    hashMap.put(nowindex, 0);
                hashMap.put(nowindex, hashMap.get(nowindex)+1);
            }
        }
        if(hashMap.keySet().size()==0)
            return hangshu;
        int ans = Integer.MAX_VALUE;
        for(Integer integer:hashMap.keySet()){
            if(hangshu - hashMap.get(integer) < ans)
                ans = hangshu - hashMap.get(integer);
        }
        return ans;
    }

    // # 690
    public int getImportance(List<Employee> employees, int id) {
        Queue<Employee> queue = new LinkedList<>();
        int ans = 0;
        HashMap<Integer, Employee> hashMap = new LinkedHashMap<>();
        for(Employee employee: employees)
            hashMap.put(employee.id, employee);
        queue.add(hashMap.get(id));
        while (!queue.isEmpty()){
            Employee employee = queue.poll();
            ans += employee.importance;
            for (Integer integer:employee.subordinates)
                queue.add(hashMap.get(integer));
        }
        return ans;
    }

    // # 740
    public int deleteAndEarn(int[] nums) {
        int Maxvalue = Arrays.stream(nums).max().getAsInt();
        int[] sum = new int[Maxvalue + 1];
        for(int i: nums)
            sum[i] += i;

        int[][] dp = new int[2][Maxvalue + 1];
        dp[0][1] = 0;
        dp[1][1] = sum[1];
        int ans = Math.max(dp[0][1], dp[1][1]);
        for(int i=2;i<Maxvalue+1;i++){
            dp[0][i] = Math.max(dp[1][i-1], dp[0][i-1]);
            if(dp[0][i] > ans)
                ans = dp[0][i];
            dp[1][i] = dp[0][i-1] + sum[i];
            if(dp[1][i] > ans)
                ans = dp[1][i];
        }
        return ans;
    }

    // # 783
    public int minDiffInBST(TreeNode root) {
        if(root==null)
            return 0;
        int ans = Integer.MAX_VALUE;
        List<Integer> list = new ArrayList<>();

        Stack<TreeNode> stack = new Stack<>();
        TreeNode node = root;
        while (node!=null || !stack.isEmpty()){
            if(node!=null){
                stack.push(node);
                node = node.left;
            }else {
                TreeNode tem = stack.pop();
                list.add(tem.val);
                node = tem.right;
            }
        }
        for(int i=1;i<list.size();i++)
            ans = Math.min(ans, list.get(i)- list.get(i-1));

        return ans;

    }

    // # 860
    public boolean lemonadeChange(int[] bills) {
        if(bills.length == 0)
            return true;
        // 第一单必须为5元账单
        if(bills[0]!=5)
            return false;
        // 定义5元，10元数量
        int five_num = 1, ten_num = 0;

        // 遍历bills
        for (int i = 1; i < bills.length; i++){
            // 5元账单，无需找零
            if(bills[i] == 5)
                five_num++;
            // 10元账单，找零一张5元
            else if(bills[i] == 10){
                if(five_num > 0) {
                    ten_num++;
                    five_num--;
                }else
                    return false;
            // 20元账单，找零一张5元，一张10元
            }else{
                if(five_num > 0 && ten_num > 0){
                    five_num --;
                    ten_num --;
                }else if(five_num >= 3){
                    five_num -= 3;
                }else
                    return false;
            }
        }

        return true;
    }

    // # 872
    public boolean leafSimilar(TreeNode root1, TreeNode root2) {
        TreeNode cur = root1;
        Stack<TreeNode> stack = new Stack<>();
        Queue<Integer> queue = new LinkedList<>();
        while (cur!=null || !stack.isEmpty()){
            while (cur!=null){
                stack.push(cur);
                cur = cur.left;
            }
            cur = stack.pop();
            if(cur.left==null && cur.right==null)
                queue.add(cur.val);
            cur = cur.right;
        }

        cur = root2;
        while (cur!=null || !stack.isEmpty()){
            while (cur!=null){
                stack.push(cur);
                cur = cur.left;
            }
            cur = stack.pop();
            if(cur.left==null && cur.right==null)
                if(queue.isEmpty() || cur.val!=queue.poll())
                    return false;
            cur = cur.right;
        }
        return queue.isEmpty()? true:false;
    }


    // # 897
    public TreeNode increasingBST(TreeNode root) {
        TreeNode ans = new TreeNode();
        TreeNode temp = ans;
        Stack<TreeNode> stack = new Stack<>();
        TreeNode cur = root;
        while (cur!=null || !stack.isEmpty()){
            while (cur!=null){
                stack.push(cur);
                cur = cur.left;
            }
            cur = stack.pop();
            temp.right = new TreeNode(cur.val);
            temp = temp.right;
            cur = cur.right;
        }
        return ans.right;

    }

    // # 938
    public int rangeSumBST(TreeNode root, int low, int high) {
        Stack<TreeNode> stack = new Stack<>();
        TreeNode cur = root;
        int ans = 0;
        while (cur!=null || !stack.isEmpty()){
            while (cur!=null){
                stack.push(cur);
                cur = cur.left;
            }
            cur = stack.pop();
            if(cur.val <= high && cur.val >= low)
                ans += cur.val;
            cur = cur.right;
        }
        return ans;
    }

    // # 1011
    public int shipWithinDays(int[] weights, int D) {
        int left=0, right = 0;
        for(int i:weights){
            if(i>left)
                left = i;
            right += i;
        }
        int ans = (left+right)/2;
        while(left<right){
            ans = (left+right)/2 ;
            System.out.println("x: " + ans);
            int count = 0;
            int sum = 0;
            for(int i:weights){
                if(sum + i <= ans){
                    sum += i;
                    continue;
                }else {
                    System.out.println("group" + count + ":" + sum);
                    sum = i;
                    count ++;
                }
            }
            count++;
            System.out.println("group" + (count-1) + ":" + sum);
            if(count > D){
                left = ans+1;
            }else {
                right = ans ;
            }
            System.out.println(left + "--" + right);

        }

        return right;
    }

    // # 1486
    public int xorOperation(int n, int start) {
        int ans = start;
        for(int i=1;i<n;i++){
            ans = ans ^ (start + 2 * i);
        }
        return ans;
    }

    // # 1720
    public int[] decode(int[] encoded, int first) {
        int[] ans = new int[encoded.length + 1];
        ans[0] = first;
        for(int i = 0; i < encoded.length; i++){
            ans[i+1] = encoded[i] ^ first;
            first = ans[i+1];
        }
        return ans;
    }

    // # 1723
    public int max_time = Integer.MAX_VALUE;
    public int minimumTimeRequired(int[] jobs, int k) {
        int[] work = new int[k];
        Arrays.sort(jobs);
        duigui(jobs, 0, k, work);
        return max_time;
    }

    public void duigui(int[] jobs, int index, int k, int[] work_k){
        if(index==jobs.length)
        {
            System.out.println(Arrays.toString(work_k));
            if(Arrays.stream(work_k).max().getAsInt() < max_time)
                max_time = Arrays.stream(work_k).max().getAsInt();
            return;
        }else {
            if(Arrays.stream(work_k).max().getAsInt() + jobs[index] > max_time)
                return;
            for(int i=0;i<work_k.length;i++){
                work_k[i] += jobs[index];
                System.out.println(index + "--" + i);
                duigui(jobs, index+1, k, work_k);
                work_k[i] -= jobs[index];
            }
            return;
        }
    }


}
