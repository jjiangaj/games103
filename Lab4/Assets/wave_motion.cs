using UnityEngine;
using System.Collections;
using System.Drawing;
using System.Numerics;
using Vector3 = UnityEngine.Vector3;

public class wave_motion : MonoBehaviour 
{
	int size 		= 100;
	float rate 		= 0.005f;
	float gamma		= 0.001f;
	float damping 	= 0.996f;

	float dt = 0.1f;
	float[,] 	old_h;
	float[,]	low_h;
	float[,]	vh;
	float[,]	b;

	bool [,]	cg_mask;
	float[,]	cg_p;
	float[,]	cg_r;
	float[,]	cg_Ap;
	bool 	tag=true;

	private int csize;

	float epsilon = 0.1f;

	Vector3[] X;

	Vector3 	cube_v = Vector3.zero;
	Vector3 	cube_w = Vector3.zero;
	UnityEngine.Matrix4x4 I_ref;							// reference inertia

	float psize = 0;

	GameObject cube;
	GameObject block;
	Collider cubeCollider;
	Collider blockCollider;
	
	// Use this for initialization
	void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.Clear ();

		X=new Vector3[size*size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			X[i*size+j].x=i*0.1f-size*0.05f;
			X[i*size+j].y=0;
			X[i*size+j].z=j*0.1f-size*0.05f;
		}

		int[] T = new int[(size - 1) * (size - 1) * 6];
		int index = 0;
		for (int i=0; i<size-1; i++)
		for (int j=0; j<size-1; j++)
		{
			T[index*6+0]=(i+0)*size+(j+0);
			T[index*6+1]=(i+0)*size+(j+1);
			T[index*6+2]=(i+1)*size+(j+1);
			T[index*6+3]=(i+0)*size+(j+0);
			T[index*6+4]=(i+1)*size+(j+1);
			T[index*6+5]=(i+1)*size+(j+0);
			index++;
		}
		mesh.vertices  = X;
		mesh.triangles = T;
		mesh.RecalculateNormals ();
		psize = mesh.bounds.size.x;

		low_h 	= new float[size,size];
		old_h 	= new float[size,size];
		vh 	  	= new float[size,size];
		b 	  	= new float[size,size];

		cg_mask	= new bool [size,size];
		cg_p 	= new float[size,size];
		cg_r 	= new float[size,size];
		cg_Ap 	= new float[size,size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			low_h[i,j]=99999;
			old_h[i,j]=0;
			vh[i,j]=0;
		}

		cube = GameObject.Find("Cube");
		cubeCollider = cube.GetComponent<Collider>();
		block = GameObject.Find("Block");
		blockCollider = block.GetComponent<Collider>();
		
		Mesh cmesh = cube.GetComponent<MeshFilter>().mesh;
		Vector3[] cvertices = mesh.vertices;

		float m=1;
		float mass=0;
		for (int i=0; i<cvertices.Length; i++) 
		{
			mass += m;
			float diag=m*cvertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*cvertices[i][0]*cvertices[i][0];
			I_ref[0, 1]-=m*cvertices[i][0]*cvertices[i][1];
			I_ref[0, 2]-=m*cvertices[i][0]*cvertices[i][2];
			I_ref[1, 0]-=m*cvertices[i][1]*cvertices[i][0];
			I_ref[1, 1]-=m*cvertices[i][1]*cvertices[i][1];
			I_ref[1, 2]-=m*cvertices[i][1]*cvertices[i][2];
			I_ref[2, 0]-=m*cvertices[i][2]*cvertices[i][0];
			I_ref[2, 1]-=m*cvertices[i][2]*cvertices[i][1];
			I_ref[2, 2]-=m*cvertices[i][2]*cvertices[i][2];
		}
		I_ref [3, 3] = 1;
		csize = cvertices.Length;
	}

	void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
	{
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			Ax[i,j]=0;
			if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
			if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
			if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
			if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
		}
	}

	float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
	{
		float ret=0;
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			ret+=x[i,j]*y[i,j];
		}
		return ret;
	}

	void Conjugate_Gradient(bool[,] mask, float[,] b, ref float[,] x, int li, int ui, int lj, int uj)
	{
		//Solve the Laplacian problem by CG.
		A_Times(mask, x, cg_r, li, ui, lj, uj);

		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
		}

		float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

		for(int k=0; k<128; k++)
		{
			if(rk_norm<1e-10f)	break;
			A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
			float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				x[i,j]   +=alpha*cg_p[i,j];
				cg_r[i,j]-=alpha*cg_Ap[i,j];
			}

			float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
			float beta=_rk_norm/rk_norm;
			rk_norm=_rk_norm;

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
			}
		}

	}
	
	UnityEngine.Quaternion QuatAdd(UnityEngine.Quaternion s, UnityEngine.Quaternion v)
	{
		UnityEngine.Quaternion q = new UnityEngine.Quaternion(s.x + v.x, s.y + v.y, s.z + v.z, s.w + v.w);
		return q;
	}

	void Shallow_Wave(ref float[,] old_h, ref float[,] h, ref float [,] new_h)
	{		
		//Step 1:
		//TODO: Compute new_h based on the shallow wave model.
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				new_h[i, j] = h[i, j] + (h[i, j] - old_h[i, j]) * damping;
				if (i - 1 > -1)   new_h[i, j] += (h[i - 1, j] - h[i, j]) * rate;
				if (i + 1 < size) new_h[i, j] += (h[i + 1, j] - h[i, j]) * rate;
				if (j - 1 > -1)   new_h[i, j] += (h[i, j - 1] - h[i, j]) * rate;
				if (j + 1 < size) new_h[i, j] += (h[i, j + 1] - h[i, j]) * rate;
			}
		}

		//Step 2: Block->Water coupling
		//TODO: for block 1, calculate low_h.
		//TODO: then set up b and cg_mask for conjugate gradient.
		//TODO: Solve the Poisson equation to obtain vh (virtual height).

		Vector3 cubeC = cube.transform.position;
		int x_min = (int)Mathf.Max(0, (cubeCollider.bounds.min.x + psize/2 + epsilon) / (psize/size));
		int x_max = (int)Mathf.Min(size, (cubeCollider.bounds.max.x + psize/2 + epsilon) / (psize/size));
		int y_min = (int)Mathf.Max(0, (cubeCollider.bounds.min.z + psize/2 + epsilon) / (psize/size));
		int y_max = (int)Mathf.Min(size, (cubeCollider.bounds.max.z + psize/2 + epsilon) / (psize/size));
		
		Vector3 [,] contact; 
		contact = new Vector3[size,size];

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (i >= x_min && i < x_max && j >= y_min && j < y_max)
				{
					Vector3 origin = new Vector3(X[i*size+j].x, -3, X[i*size + j].z);
					Ray ray = new Ray(origin, transform.up);
					RaycastHit hit;
					if (cubeCollider.Raycast(ray, out hit, 10))
					{
						if (hit.point.y < h[i, j])
						{
							low_h[i, j] = hit.point.y;
							b[i, j] = (new_h[i, j] - low_h[i, j]) * 200;
							cg_mask[i, j] = true;

							//Vector3 l = hit.point - cubeC;
							//l /= l.magnitude;            // unit vector from center to contact point
							contact[i, j] = hit.point;
						}
					}
					else
					{
						vh[i, j] = 0;
						cg_mask[i, j] = false;
					}
				}
				else
				{
					vh[i, j] = 0;
					cg_mask[i, j] = false;
				}
			}
		}
		Conjugate_Gradient(cg_mask, b,  ref vh, x_min, x_max, y_min, y_max);
		
		Vector3 Fg = new Vector3(0, -2.6f, 0);
		Vector3 torque = new Vector3();
		Vector3 force = new Vector3();
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (contact[i, j] != Vector3.zero)
				{
					Vector3 le = contact[i,j] - cubeC;
					le /= le.magnitude;

					Vector3 fi = Vector3.Dot  (le, vh[i, j] * transform.up) * le;
					force += fi;

					Vector3 ti = Vector3.Cross(le, vh[i, j] * transform.up);
					torque += ti;
				}
			}
		}

		UnityEngine.Quaternion q = cube.transform.rotation;
		Vector3 w = I_ref.inverse.MultiplyPoint(torque);
		UnityEngine.Quaternion dq = new UnityEngine.Quaternion(dt*w.x/2, dt*w.y/2, dt*w.z/2, 0.0f);
		q = QuatAdd(q, dq * q);
		q = UnityEngine.Quaternion.Normalize(q);
		
		Vector3 pos = cubeC + dt * dt * (Fg + force/(csize)) / 2;
		cube.transform.position = pos;
		cube.transform.rotation = q;

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (i - 1 > -1)   new_h[i, j] += (vh[i - 1, j] - vh[i, j]) * rate * gamma;
				if (i + 1 < size) new_h[i, j] += (vh[i + 1, j] - vh[i, j]) * rate * gamma;
				if (j - 1 > -1)   new_h[i, j] += (vh[i, j - 1] - vh[i, j]) * rate * gamma;
				if (j + 1 < size) new_h[i, j] += (vh[i, j + 1] - vh[i, j]) * rate * gamma;
			}
		}
		
		
		//TODO: for block 2, calculate low_h.
		//TODO: then set up b and cg_mask for conjugate gradient.
		//TODO: Solve the Poisson equation to obtain vh (virtual height).
		Vector3 blockC = block.transform.position;
		int x2_min = (int)Mathf.Max(0, (blockCollider.bounds.min.x + psize/2 + epsilon) / (psize/size));
		int x2_max = (int)Mathf.Min(size, (blockCollider.bounds.max.x + psize/2 + epsilon) / (psize/size));
		int y2_min = (int)Mathf.Max(0, (blockCollider.bounds.min.z + psize/2 + epsilon) / (psize/size));
		int y2_max = (int)Mathf.Min(size, (blockCollider.bounds.max.z + psize/2 + epsilon) / (psize/size));

		Vector3 [,] contact2; 
		contact2 = new Vector3[size,size];

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (i >= x2_min && i < x2_max && j >= y2_min && j < y2_max)
				{
					Vector3 origin = new Vector3(X[i*size+j].x, -3, X[i*size + j].z);
					Ray ray = new Ray(origin, transform.up);
					RaycastHit hit;
					if (blockCollider.Raycast(ray, out hit, 10))
					{
						if (hit.point.y < h[i, j])
						{
							low_h[i, j] = hit.point.y;
							b[i, j] = (new_h[i, j] - low_h[i, j]) * 200;
							cg_mask[i, j] = true;
							contact2[i, j] = hit.point;
						}
					}
					else
					{
						vh[i, j] = 0;
						cg_mask[i, j] = false;
					}
				}
				else
				{
					vh[i, j] = 0;
					cg_mask[i, j] = false;
				}
			}
		}
		Conjugate_Gradient(cg_mask, b,  ref vh, x2_min, x2_max, y2_min, y2_max);
		
		Vector3 torque2 = new Vector3();
		Vector3 force2 = new Vector3();
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (contact2[i, j] != Vector3.zero)
				{
					Vector3 le = contact2[i,j] - blockC;
					le /= le.magnitude;

					Vector3 fi = Vector3.Dot  (le, vh[i, j] * transform.up) * le;
					force2 += fi;

					Vector3 ti = Vector3.Cross(le, vh[i, j] * transform.up);
					torque2 += ti;
				}
			}
		}

		UnityEngine.Quaternion q2 = block.transform.rotation;
		Vector3 w2 = I_ref.inverse.MultiplyPoint(torque2);
		UnityEngine.Quaternion dq2 = new UnityEngine.Quaternion(dt*w2.x/2, dt*w2.y/2, dt*w2.z/2, 0.0f);
		q2 = QuatAdd(q2, dq2 * q2);
		q2 = UnityEngine.Quaternion.Normalize(q2);
		
		Vector3 pos2 = blockC + dt * dt * (Fg + force2/(csize)) / 2;
		block.transform.position = pos2;
		block.transform.rotation = q2;
		
		//TODO: Diminish vh.
		//TODO: Update new_h by vh.

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (i - 1 > -1)   new_h[i, j] += (vh[i - 1, j] - vh[i, j]) * rate * gamma;
				if (i + 1 < size) new_h[i, j] += (vh[i + 1, j] - vh[i, j]) * rate * gamma;
				if (j - 1 > -1)   new_h[i, j] += (vh[i, j - 1] - vh[i, j]) * rate * gamma;
				if (j + 1 < size) new_h[i, j] += (vh[i, j + 1] - vh[i, j]) * rate * gamma;
			}
		}

		//Step 3
		//TODO: old_h <- h; h <- new_h;
		old_h = h.Clone() as float[,];
		h = new_h.Clone() as float[,];

		//Step 4: Water->Block coupling.
		//More TODO here.
		
	}
	

	// Update is called once per frame
	void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		X    = mesh.vertices;
		float[,] new_h = new float[size, size];
		float[,] h     = new float[size, size];

		//TODO: Load X.y into h.
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				h[i, j] = X[i * size + j].y;
			}
		}

		if (Input.GetKeyDown ("r")) 
		{
			//TODO: Add random water.
			int xwater = Random.Range(1, size-1);
			int ywater = Random.Range(1, size-1);
			float r = Random.Range(0.1f, 0.2f);
			
			
			h[xwater, ywater] += r;
			h[xwater, ywater - 1] -= r / 8;
			h[xwater, ywater + 1] -= r / 8;
			h[xwater - 1, ywater] -= r / 8;
			h[xwater + 1, ywater] -= r / 8;
			
			h[xwater - 1, ywater + 1] -= r / 8;
			h[xwater - 1, ywater - 1] -= r / 8;
			h[xwater + 1, ywater + 1] -= r / 8;
			h[xwater + 1, ywater - 1] -= r / 8;
		}
	
		for(int l=0; l<8; l++)
		{
			Shallow_Wave(ref old_h, ref h, ref new_h);
		}

		//TODO: Store h back into X.y and recalculate normal.
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				X[i * size + j].y = h[i,j];
			}
		}
		
		mesh.vertices  = X;
		mesh.RecalculateNormals ();

	}
}
