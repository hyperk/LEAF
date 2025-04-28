#include "LeafUtility.hh"

bool ContainsTrueVtx(std::vector<std::vector<double>>* tRecoVtxPos, std::vector<double> fTrueVtxPos)
{
	for (unsigned int i = 0; i < tRecoVtxPos->size(); i++)
	{
		double distance = Distance3D((*tRecoVtxPos)[i], fTrueVtxPos);
		if (distance <= fSearchVtxStep)
		{
			return true;
		}
	}
	return false;
}

//* to know if the predicted vertex and the true vertex belong to the same candidate
//? doesn't work because the radius overlap between candidates
bool CorrectCandidate(std::vector<std::vector<double>>* candidates, std::vector<double> point, std::vector<double> fTrueVtxPos)
{
	for(unsigned int i = 0; i < candidates->size(); i++)
	{
		if (Distance3D((*candidates)[i], fTrueVtxPos) <= fSearchVtxStep && Distance3D((*candidates)[i], point) <= fSearchVtxStep)
		{
			return true;
		}
	}
	return false;
}

double Distance3D(std::vector<double> point1, std::vector<double> point2)
{
	return sqrt(pow(point1[0] - point2[0], 2) + 
				pow(point1[1] - point2[1], 2) + 
				pow(point1[2] - point2[2], 2));
}