using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class PBD_model : MonoBehaviour {

    int[] E;
    float[] L;
    Vector3[] V;

    const float t = 0.0333f;
    const float t_neg = 1 / t;
    const float t2_neg = 1 / (t * t);
    const float mass = 1;
    const float damping = 0.99f;
    const float rho = 0.995f;
    const float spring_k = 8000;
    const float spring_k4 = spring_k * 4;
    static readonly HashSet<int> fixedPoint = new HashSet<int> { 0, 20 };
    const int N = 21;       // 将 mesh 重构为 20*20 的网格



    #region Initialization
    // Use this for initialization
    void Start() {
        // 隐藏之前的 canvas 文字
        GameObject.Find("Canvas").GetComponent<Canvas>().enabled = false;

        Mesh mesh = GetComponent<MeshFilter>().mesh;

        // Resize the mesh.
        int n = 21;
        Vector3[] X = new Vector3[n * n];
        Vector2[] UV = new Vector2[n * n];
        int[] T = new int[(n - 1) * (n - 1) * 6];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                X[j * n + i] = new Vector3(5 - 10.0f * i / (n - 1), 0, 5 - 10.0f * j / (n - 1));
                UV[j * n + i] = new Vector3(i / (n - 1.0f), j / (n - 1.0f));
            }
        int t = 0;
        for (int j = 0; j < n - 1; j++)
            for (int i = 0; i < n - 1; i++) {
                T[t * 6 + 0] = j * n + i;
                T[t * 6 + 1] = j * n + i + 1;
                T[t * 6 + 2] = (j + 1) * n + i + 1;
                T[t * 6 + 3] = j * n + i;
                T[t * 6 + 4] = (j + 1) * n + i + 1;
                T[t * 6 + 5] = (j + 1) * n + i;
                t++;
            }
        mesh.vertices = X;
        mesh.triangles = T;
        mesh.uv = UV;
        mesh.RecalculateNormals();

        // Construct the original edge list
        int[] _E = new int[T.Length * 2];
        for (int i = 0; i < T.Length; i += 3) {
            _E[i * 2 + 0] = T[i + 0];
            _E[i * 2 + 1] = T[i + 1];
            _E[i * 2 + 2] = T[i + 1];
            _E[i * 2 + 3] = T[i + 2];
            _E[i * 2 + 4] = T[i + 2];
            _E[i * 2 + 5] = T[i + 0];
        }
        // Reorder the original edge list
        for (int i = 0; i < _E.Length; i += 2)
            if (_E[i] > _E[i + 1])
                Swap(ref _E[i], ref _E[i + 1]);
        // Sort the original edge list using quicksort
        Quick_Sort(ref _E, 0, _E.Length / 2 - 1);

        int e_number = 0;
        for (int i = 0; i < _E.Length; i += 2)
            if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1])
                e_number++;

        E = new int[e_number * 2];
        for (int i = 0, e = 0; i < _E.Length; i += 2)
            if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1]) {
                E[e * 2 + 0] = _E[i + 0];
                E[e * 2 + 1] = _E[i + 1];
                e++;
            }

        L = new float[E.Length / 2];
        for (int e = 0; e < E.Length / 2; e++) {
            int i = E[e * 2 + 0];
            int j = E[e * 2 + 1];
            L[e] = (X[i] - X[j]).magnitude;
        }

        V = new Vector3[X.Length];
        for (int i = 0; i < X.Length; i++)
            V[i] = new Vector3(0, 0, 0);
    }

    void Quick_Sort(ref int[] a, int l, int r) {
        int j;
        if (l < r) {
            j = Quick_Sort_Partition(ref a, l, r);
            Quick_Sort(ref a, l, j - 1);
            Quick_Sort(ref a, j + 1, r);
        }
    }

    int Quick_Sort_Partition(ref int[] a, int l, int r) {
        int pivot_0, pivot_1, i, j;
        pivot_0 = a[l * 2 + 0];
        pivot_1 = a[l * 2 + 1];
        i = l;
        j = r + 1;
        while (true) {
            do ++i; while (i <= r && (a[i * 2] < pivot_0 || a[i * 2] == pivot_0 && a[i * 2 + 1] <= pivot_1));
            do --j; while (a[j * 2] > pivot_0 || a[j * 2] == pivot_0 && a[j * 2 + 1] > pivot_1);
            if (i >= j) break;
            Swap(ref a[i * 2], ref a[j * 2]);
            Swap(ref a[i * 2 + 1], ref a[j * 2 + 1]);
        }
        Swap(ref a[l * 2 + 0], ref a[j * 2 + 0]);
        Swap(ref a[l * 2 + 1], ref a[j * 2 + 1]);
        return j;
    }

    void Swap(ref int a, ref int b) {
        int temp = a;
        a = b;
        b = temp;
    }
    #endregion Initialization

    void Strain_Limiting() {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;

        // [2.b] Apply PBD here.

        // 初始化
        int length = X.Length;
        Vector3[] sum_X = new Vector3[length];
        int[] sum_n = new int[length];
        for (int i = 0; i < length; ++i) {
            sum_X[i] = Vector3.zero;
            sum_n[i] = 0;
        }

        // 限制
        length = E.Length / 2;
        for (int i = 0; i < length; ++i) {
            int tmp = 2 * i;
            int l = E[tmp], r = E[tmp + 1];
            Vector3 dis = (X[l] - X[r]).normalized * L[i];
            Vector3 add = X[l] + X[r];
            sum_X[l] += 0.5F * (add + dis);
            sum_X[r] += 0.5F * (add - dis);
            ++sum_n[l];
            ++sum_n[r];
        }

        // 更新
        length = X.Length;
        for (int i = 0; i < length; ++i) {
            if (fixedPoint.Contains(i)) { continue; }

            Vector3 new_X = (0.2F * X[i] + sum_X[i]) / (0.2F + sum_n[i]);
            V[i] += t_neg * (new_X - X[i]);
            X[i] = new_X;
        }

        mesh.vertices = X;
    }

    void Collision_Handling() {
        // [1.e] Sphere Collision
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;
        GameObject sphere = GameObject.Find("Sphere");
        Vector3 center = sphere.transform.position;
        const float r = 2.7F;
        const float r2 = r * r;

        // Handle collision
        for (int i = 0; i < X.Length; ++i) {
            if (fixedPoint.Contains(i)) { continue; }
            Vector3 disV = X[i] - center;
            float dis2 = disV.sqrMagnitude;
            if (dis2 < r2) {
                // collision
                float dis = Mathf.Sqrt(dis2);
                Vector3 disR = r * disV / dis;
                V[i] += t_neg * (center + disR - X[i]);
                X[i] = center + disR;
            }
        }

        mesh.vertices = X;
    }


    // Update is called once per frame
    void Update() {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;

        // [2.a]
        int length = X.Length;


        Vector3 g = new Vector3(0, -9.8F, 0);

        // 速度阻尼, 更新速度, 更新位置
        for (int i = 0; i < length; ++i) {
            if (fixedPoint.Contains(i)) { continue; }
            Vector3 v = V[i] * damping + g * t;
            V[i] *= damping;
            V[i] = v;
            X[i] += v * t;
        }

        mesh.vertices = X;

        for (int l = 0; l < 32; l++) {
            Strain_Limiting();
        }

        Collision_Handling();

        mesh.RecalculateNormals();

    }
}